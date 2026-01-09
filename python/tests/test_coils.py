"""Tests for coils visualization module."""

from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pytest

from libneo.coils import CoilsFile


@pytest.fixture
def simple_coils_file():
    """Create a temporary coils file for testing."""
    content = """periods 5
mirror NIL
begin filament
0.0 1.0 0.0 1000.0
0.1 1.0 0.1 1000.0
0.2 1.0 0.0 1000.0
end
begin filament
0.0 -1.0 0.0 -1000.0
0.1 -1.0 0.1 -1000.0
0.2 -1.0 0.0 -1000.0
end
"""
    with NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        return Path(f.name)


def test_coilsfile_from_file(simple_coils_file):
    """Test CoilsFile.from_file from libneo."""
    coils = CoilsFile.from_file(str(simple_coils_file))

    assert len(coils.filaments) == 2
    assert coils.periods == 5
    assert coils.mirror == "NIL"

    fil1 = coils.filaments[0]
    assert fil1.coords.shape == (3, 3)
    assert np.allclose(fil1.current, 1000.0)


def test_coil_to_tube_mesh(simple_coils_file):
    """Test creating tube mesh from coil coords."""
    pytest.importorskip("trimesh")

    from gecko.coils import coil_to_tube_mesh

    coils = CoilsFile.from_file(str(simple_coils_file))
    fil = coils.filaments[0]

    mesh = coil_to_tube_mesh(fil.coords, width=0.1, height=0.2)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_coils_to_mesh(simple_coils_file):
    """Test converting all coils to combined mesh."""
    pytest.importorskip("trimesh")

    from gecko.coils import coils_to_mesh

    coils = CoilsFile.from_file(str(simple_coils_file))
    mesh = coils_to_mesh(coils, width=0.1, height=0.2)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
