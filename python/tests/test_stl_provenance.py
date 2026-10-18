"""Exported STL provenance must identify the producing source revision."""

import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest
import trimesh

from gecko.cli import main
from gecko.provenance import geometry_version
from gecko.stl import convert_stl_units, export_mesh_to_stl


def read_provenance(path):
    data = json.loads(path.with_suffix(path.suffix + ".json").read_text())
    assert data["stl_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert data["gecko_version"]
    return data


def test_export_revision_and_geometry(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[2]
    expected = subprocess.check_output(
        ["git", "describe", "--tags", "--always", "--dirty"], cwd=root, text=True
    ).strip()
    monkeypatch.chdir(tmp_path)  # caller need not be inside the source repository
    mesh = trimesh.creation.box(extents=(1, 2, 3))
    path = export_mesh_to_stl(mesh, tmp_path / "nested" / "box.stl")
    assert read_provenance(path)["git_describe"] == expected
    np.testing.assert_allclose(trimesh.load_mesh(path).extents, mesh.extents)


def test_conversion_updates_provenance(tmp_path):
    source = export_mesh_to_stl(
        trimesh.creation.box(extents=(1000, 2000, 3000)), tmp_path / "input.stl"
    )
    output = convert_stl_units(
        source, tmp_path / "output.stl", from_units="mm", to_units="m"
    )
    before, after = read_provenance(source), read_provenance(output)
    assert before["stl_sha256"] != after["stl_sha256"]
    assert before["git_describe"] == after["git_describe"]
    np.testing.assert_allclose(trimesh.load_mesh(output).extents, [1, 2, 3])


@pytest.mark.parametrize("offset", ["0", "0.05"])
def test_wall_cli_stamps_stl(tmp_path, offset):
    wout = Path(__file__).parent / "data" / "alpes" / "wout.nc"
    output = tmp_path / "wall.stl"
    assert (
        main(
            [
                "wall-from-vmec",
                str(wout),
                "--offset",
                offset,
                "--n-theta",
                "32",
                "--n-zeta",
                "24",
                "-o",
                str(output),
            ]
        )
        == 0
    )
    assert read_provenance(output)["git_describe"] == geometry_version()["git_describe"]


def test_git_unavailable(monkeypatch):
    import sys

    monkeypatch.setitem(sys.modules, "gecko._build_revision", None)

    def missing_git(*args, **kwargs):
        raise FileNotFoundError("git unavailable")

    monkeypatch.setattr(subprocess, "check_output", missing_git)
    result = geometry_version()
    assert result["gecko_version"]
    assert result["git_describe"] is None


def test_packaged_revision_without_git(monkeypatch):
    import sys
    from types import SimpleNamespace

    monkeypatch.setitem(
        sys.modules,
        "gecko._build_revision",
        SimpleNamespace(GIT_DESCRIBE="v0.1.0-2-g1234567"),
    )

    def missing_git(*args, **kwargs):
        raise FileNotFoundError("git unavailable")

    monkeypatch.setattr(subprocess, "check_output", missing_git)
    assert geometry_version()["git_describe"] == "v0.1.0-2-g1234567"
