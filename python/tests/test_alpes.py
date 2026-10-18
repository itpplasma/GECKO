"""Behavioral tests using the open ALPES geometry (CC BY 4.0)."""

from pathlib import Path

import numpy as np
import pytest
from libneo.coils import CoilsFile, Filament
from libneo.vmec import VMECGeometry
from shapely.geometry import Polygon

from gecko.coils import coil_to_tube_mesh, coils_to_mesh
from gecko.collision import detect_coil_collisions
from gecko.wall import offset_wall_mesh

DATA = Path(__file__).parent / "data" / "alpes"


@pytest.fixture(scope="module")
def alpes_coils():
    filaments = []
    for i in range(12):
        coords = np.loadtxt(DATA / f"coil_coordinates{i}.txt")
        filaments.append(Filament(coords=coords, current=np.ones(len(coords))))
    return CoilsFile(periods=2, mirror="NIL", filaments=filaments)


@pytest.fixture(scope="module")
def offset_wall():
    return offset_wall_mesh(DATA / "wout.nc", 0.05, n_theta=128, n_zeta=96)


def test_offset_wall_closed_and_outward(offset_wall):
    assert np.isfinite(offset_wall.vertices).all()
    assert offset_wall.is_watertight
    assert offset_wall.is_winding_consistent
    assert offset_wall.volume > 0
    assert offset_wall.euler_number == 0  # a closed toroidal wall
    assert np.all(offset_wall.area_faces > 0)


def test_offset_wall_lcfs_clearance(offset_wall):
    # The offset is defined in RZ planes, not along the 3D surface normal.
    # Evaluate the LCFS independently at 8x the wall's poloidal resolution.
    geom = VMECGeometry.from_file(str(DATA / "wout.nc"))
    theta = np.linspace(0, 2 * np.pi, 1024, endpoint=False)
    rings = offset_wall.vertices.reshape(96, 128, 3)
    for i, zeta in enumerate(np.linspace(0, 2 * np.pi, 96, endpoint=False)):
        r, z, _ = geom.coords(-1, theta, zeta)
        lcfs = Polygon(np.column_stack((r, z)))
        ring = rings[i]
        wall = Polygon(np.column_stack((np.hypot(ring[:, 0], ring[:, 1]), ring[:, 2])))
        assert wall.is_valid
        assert wall.contains(lcfs)
        # Allow 2 mm for the polygonal buffer and resampling discretization.
        assert 0.048 < wall.boundary.distance(lcfs.boundary) < 0.052


def test_coil_tubes(alpes_coils):
    for filament in alpes_coils.filaments:
        mesh = coil_to_tube_mesh(filament.coords, width=0.01, height=0.02)
        assert np.isfinite(mesh.vertices).all()
        assert mesh.is_watertight
        assert mesh.is_winding_consistent
        assert np.all(mesh.area_faces > 0)
        # Every rectangular profile has the requested size around its centerline.
        profiles = mesh.vertices.reshape(-1, 4, 3)
        np.testing.assert_allclose(profiles.mean(axis=1), filament.coords, atol=1e-14)
        edges = np.linalg.norm(np.roll(profiles, -1, axis=1) - profiles, axis=2)
        np.testing.assert_allclose(
            edges, np.tile([0.01, 0.02, 0.01, 0.02], (len(filament.coords), 1))
        )
    combined = coils_to_mesh(alpes_coils, width=0.01, height=0.02)
    assert len(combined.split()) == 12


def test_coils_do_not_collide(alpes_coils):
    result = detect_coil_collisions(alpes_coils, width=0.01, height=0.02)
    assert result.n_coils == 12
    assert result.collision_pairs == []
    assert not result.has_collisions


def test_coils_detect_overlap(alpes_coils):
    # Overlap two copies of a real ALPES filament with a small translation.
    # This checks the contact path, not just bounding-box rejection.
    first = alpes_coils.filaments[0]
    shifted = Filament(first.coords + [0.001, 0, 0], first.current.copy())
    pair = CoilsFile(periods=2, mirror="NIL", filaments=[first, shifted])
    result = detect_coil_collisions(pair, width=0.01, height=0.02)
    assert result.n_coils == 2
    assert result.has_collisions
    assert result.collision_pairs == [(0, 1)]
    contacts = result.collisions[0].contact_points
    assert contacts.shape[0] > 0
    assert contacts.shape[1] == 3
    assert np.isfinite(contacts).all()
