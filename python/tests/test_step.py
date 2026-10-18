"""STEP import behavior and open-data geometry/field round trip."""

import os
from pathlib import Path

import numpy as np
import pytest

cq = pytest.importorskip("cadquery")

from gecko.cli import main
from gecko.step import (
    StepCoilError,
    _section_centers,
    coils_from_step,
    extract_step_centerlines,
    write_step_coils,
)
from libneo.coils import CoilsFile

from .step_helpers import field, loft_coil, nearest_on_polyline


@pytest.fixture(scope="module")
def alpes_step(tmp_path_factory):
    original = (
        np.loadtxt(Path(__file__).parent / "data/alpes/coil_coordinates0.txt") * 330.0
    )
    shape = loft_coil(original)
    path = tmp_path_factory.mktemp("step") / "alpes.step"
    cq.exporters.export(shape, str(path))
    return path, original * 1e-3


@pytest.fixture(scope="module")
def ring_step(tmp_path_factory):
    angle = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    points = np.column_stack(
        [100 * np.cos(angle), 8 * np.sin(2 * angle), 100 * np.sin(angle)]
    )
    shape = loft_coil(points, width=8, height=12)
    directory = tmp_path_factory.mktemp("ring")
    path = directory / "rings.step"
    assembly = cq.Compound.makeCompound([shape, shape.translate((300, 0, 0))])
    cq.exporters.export(assembly, str(path))
    return path, shape


def test_alpes_roundtrip(alpes_step, tmp_path):
    path, original = alpes_step
    output = tmp_path / "coils.alpes"
    result = write_step_coils(path, output, currents=1200, n_points=1024)
    extracted = result.filaments[0].coords
    assert extracted.shape == (1025, 3)
    np.testing.assert_array_equal(extracted[0], extracted[-1])
    # Isolate CAD section-center extraction before uniform arc-length resampling.
    solid = cq.importers.importStep(str(path)).val().Solids()[0]
    centers = _section_centers(solid) * 1e-3
    assert centers.shape == original.shape
    center_distances = np.linalg.norm(centers[:, None, :] - original, axis=2)
    extraction_error = max(
        center_distances.min(axis=0).max(), center_distances.min(axis=1).max()
    )
    assert extraction_error < 1e-9
    # Resampled vertices lie on source segments, but new chords cut corners.
    # Sample the reverse distance on 32 positions per source segment.
    dense = np.concatenate(
        [
            original + f * (np.roll(original, -1, axis=0) - original)
            for f in np.linspace(0, 1, 32, endpoint=False)
        ]
    )
    projected, resampled_to_source = nearest_on_polyline(extracted[:-1], original)
    _, source_to_resampled = nearest_on_polyline(dense, extracted[:-1])
    assert resampled_to_source.max() < 1e-9
    assert source_to_resampled.max() < 2e-4  # 0.2 mm chord error
    assert np.all(np.linalg.norm(np.diff(extracted, axis=0), axis=1) > 0)
    assert np.linalg.norm(extracted[0] - extracted[-2]) < 0.01
    # CAD has no physical current sense; compare with the source orientation.
    source_area = np.cross(original, np.roll(original, -1, axis=0)).sum(axis=0)
    result_area = np.cross(extracted[:-1], extracted[1:]).sum(axis=0)
    if np.dot(source_area, result_area) < 0:
        original = original[::-1]
    center_area = np.cross(centers, np.roll(centers, -1, axis=0)).sum(axis=0)
    if np.dot(center_area, result_area) < 0:
        centers = centers[::-1]
    grid = np.array(np.meshgrid(*[[-0.08, 0.0, 0.08]] * 3)).reshape(3, -1).T
    _, clearance = nearest_on_polyline(grid, original)
    assert clearance.min() > 0.03
    source_field = field(grid, original, 1200)
    center_field = field(grid, centers, 1200)
    result_field = field(grid, extracted[:-1], 1200)
    extraction_relative = np.linalg.norm(center_field - source_field, axis=1) / (
        np.linalg.norm(source_field, axis=1)
    )
    resampling_relative = np.linalg.norm(result_field - center_field, axis=1) / (
        np.linalg.norm(center_field, axis=1)
    )
    assert extraction_relative.max() < 1e-9
    assert resampling_relative.max() < 2e-3
    relative = np.linalg.norm(result_field - source_field, axis=1) / np.linalg.norm(
        source_field, axis=1
    )
    assert relative.max() < 2e-3  # 0.2% at all 27 sample points
    loaded = CoilsFile.from_file(str(output))
    assert len(loaded.filaments) == 1
    # libneo strips the zero-current closing marker. Closure is explicit in file.
    np.testing.assert_allclose(loaded.filaments[0].coords, extracted[:-1], atol=1e-14)
    np.testing.assert_allclose(loaded.filaments[0].current, 1200)
    lines = output.read_text().splitlines()
    closing = lines[-2].split()
    np.testing.assert_allclose(np.array(closing[:3], float), extracted[0])
    assert closing[3:] == ["0", "1", "coil_1"]
    artifacts = Path(os.environ.get("GECKO_STEP_ARTIFACT_DIR", tmp_path))
    artifacts.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        artifacts / "centerline.csv",
        np.column_stack([projected, extracted[:-1], resampled_to_source]),
        delimiter=",",
        header="source_projection_x_m,source_projection_y_m,source_projection_z_m,"
        "resampled_x_m,resampled_y_m,resampled_z_m,resampled_to_source_distance_m",
        comments="",
    )
    np.savetxt(
        artifacts / "field_errors.csv",
        np.column_stack(
            [
                grid,
                source_field,
                center_field,
                result_field,
                extraction_relative,
                resampling_relative,
                relative,
            ]
        ),
        delimiter=",",
        header="x_m,y_m,z_m,source_Bx_T,source_By_T,source_Bz_T,"
        "section_center_Bx_T,section_center_By_T,section_center_Bz_T,"
        "resampled_Bx_T,resampled_By_T,resampled_Bz_T,"
        "extraction_relative_field_error,resampling_relative_field_error,"
        "total_relative_field_error",
        comments="",
    )
    np.savetxt(
        artifacts / "statistics.csv",
        [
            [
                extraction_error,
                resampled_to_source.max(),
                source_to_resampled.max(),
                extraction_relative.max(),
                resampling_relative.max(),
                np.sqrt(np.mean(resampling_relative**2)),
                relative.max(),
                1e-9,
                2e-4,
                1e-9,
                2e-3,
            ]
        ],
        delimiter=",",
        header="max_extraction_section_center_error_m,"
        "max_resampled_to_source_distance_m,max_source_to_resampled_distance_m,"
        "max_extraction_relative_field_error,max_resampling_relative_field_error,"
        "rms_resampling_relative_field_error,max_total_relative_field_error,"
        "extraction_tolerance_m,resampling_chord_tolerance_m,"
        "extraction_field_tolerance,resampling_and_total_field_tolerance",
        comments="",
    )


def test_units_currents_orientation_and_multiple_solids(ring_step, tmp_path):
    path, _ = ring_step
    a = coils_from_step(path, currents=[100, -200], n_points=96)
    b = coils_from_step(path, currents=10, units="m", n_points=96, reverse=[1])
    assert len(a.filaments) == 2
    for i in range(2):
        assert a.filaments[i].coords.shape == (97, 3)
        np.testing.assert_allclose(a.filaments[i].current, [100, -200][i])
        expected = b.filaments[i].coords[:: (-1 if i == 1 else 1)]
        np.testing.assert_allclose(a.filaments[i].coords, expected, atol=1e-10)
    assert a.filaments[1].coords[:, 0].mean() > 0.2
    assert (
        main(
            [
                "coils-from-step",
                str(path),
                "-o",
                str(tmp_path / "coils.test"),
                "--currents",
                "100",
                "-200",
                "--n-points",
                "96",
            ]
        )
        == 0
    )
    loaded = CoilsFile.from_file(str(tmp_path / "coils.test"))
    assert len(loaded.filaments) == 2
    np.testing.assert_allclose(loaded.filaments[1].current, -200)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"currents": [1, 2, 3]}, "one finite current"),
        ({"currents": float("nan")}, "one finite current"),
        ({"currents": 1, "reverse": [2]}, "zero-based"),
        ({"currents": 1, "n_points": 3}, "integer >= 8"),
        ({"currents": 1, "n_points": 12.5}, "integer >= 8"),
        ({"currents": 1, "units": "cm"}, "units must"),
    ],
)
def test_bad_options(ring_step, kwargs, match):
    with pytest.raises(ValueError, match=match):
        coils_from_step(ring_step[0], **kwargs)


def test_missing_and_corrupt_files(tmp_path, capsys):
    with pytest.raises(FileNotFoundError):
        extract_step_centerlines(tmp_path / "absent.step")
    path = tmp_path / "bad.step"
    path.write_text("not a STEP file")
    with pytest.raises(StepCoilError, match="Cannot read STEP"):
        extract_step_centerlines(path)
    assert (
        main(
            [
                "coils-from-step",
                str(path),
                "-o",
                str(tmp_path / "coils.bad"),
                "--currents",
                "1",
            ]
        )
        == 1
    )
    assert "STEP import failed" in capsys.readouterr().err
    assert not (tmp_path / "coils.bad").exists()


@pytest.mark.parametrize("kind", ["open", "capped", "unsupported"])
def test_reject_non_winding_packs(kind, tmp_path, ring_step):
    if kind == "open":
        shape = cq.Shell.makeShell(ring_step[1].Faces()[:-1])
        match = "closed solids"
    elif kind == "capped":
        shape = cq.Workplane().box(10, 10, 10).val()
        match = "closed loft grid"
    else:
        shape = cq.Solid.makeTorus(100, 10)
        match = "Unsupported"
    path = tmp_path / f"{kind}.step"
    cq.exporters.export(shape, str(path))
    with pytest.raises(StepCoilError, match=match):
        extract_step_centerlines(path)


def test_reject_smooth_bspline_sweep(tmp_path):
    # A valid closed rectangular sweep has no segmented section stations.
    angle = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    points = [cq.Vector(100 * np.cos(t), 0, 100 * np.sin(t)) for t in angle]
    edge = cq.Edge.makeSpline(points, periodic=True)
    rail = cq.Wire.assembleEdges([edge])
    plane = cq.Plane(origin=edge.startPoint(), normal=edge.tangentAt(0))
    shape = cq.Workplane(plane).rect(8, 12).sweep(rail, isFrenet=True).val()
    assert shape.isValid()
    assert all(shell.Closed() for shell in shape.Shells())
    assert any(face.geomType() == "BSPLINE" for face in shape.Faces())
    path = tmp_path / "smooth_sweep.step"
    cq.exporters.export(shape, str(path))
    with pytest.raises(StepCoilError, match="Unsupported topology"):
        extract_step_centerlines(path)


def test_zero_current_coil_omitted(ring_step, tmp_path):
    output = tmp_path / "coils.zero"
    coils = write_step_coils(ring_step[0], output, currents=[0, -200], n_points=32)
    assert len(coils.filaments) == 2
    loaded = CoilsFile.from_file(str(output))
    assert len(loaded.filaments) == 1
    np.testing.assert_allclose(loaded.filaments[0].current, -200)


def test_field_oracle_against_square_loop():
    # Four finite straight sides, evaluated at the center of a 2 m square.
    square = np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]])
    np.testing.assert_allclose(
        field(np.zeros((1, 3)), square, 1), [[0, 0, 4 * np.sqrt(2) * 1e-7]]
    )


def test_reject_self_intersecting_solid(tmp_path):
    # A sewn pentagram tube is topologically closed, but crosses itself.
    angle = np.arange(5) * 4 * np.pi / 5
    points = np.column_stack([100 * np.cos(angle), np.zeros(5), 100 * np.sin(angle)])
    shape = loft_coil(points, width=8, height=12, fuse=False)
    assert all(shell.Closed() for shell in shape.Shells())
    path = tmp_path / "self_intersecting.step"
    cq.exporters.export(shape, str(path))
    with pytest.raises(StepCoilError, match="[Ss]elf-intersecting"):
        extract_step_centerlines(path)


def test_all_zero_output_rejected(ring_step, tmp_path):
    output = tmp_path / "coils.zero"
    with pytest.raises(ValueError, match="At least one nonzero"):
        write_step_coils(ring_step[0], output, currents=0)
    assert not output.exists()
