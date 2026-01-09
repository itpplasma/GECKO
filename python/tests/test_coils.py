"""Tests for coils module."""

from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pytest

from gecko.coils import CoilsFile, Filament, parse_coils_file


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


def test_parse_coils_file(simple_coils_file):
    """Test parsing a simple coils file."""
    coils = parse_coils_file(simple_coils_file)

    assert coils.periods == 5
    assert coils.mirror == "NIL"
    assert len(coils.filaments) == 2

    fil1 = coils.filaments[0]
    assert fil1.coords.shape == (3, 3)
    assert np.allclose(fil1.current, 1000.0)

    fil2 = coils.filaments[1]
    assert np.allclose(fil2.current, -1000.0)


def test_coilsfile_from_file(simple_coils_file):
    """Test CoilsFile.from_file."""
    coils = CoilsFile.from_file(simple_coils_file)

    assert coils.n_coils == 2
    assert coils.total_points == 6


def test_coilsfile_write_roundtrip(simple_coils_file, tmp_path):
    """Test that write produces valid output."""
    coils = CoilsFile.from_file(simple_coils_file)
    output = tmp_path / "coils_out.txt"
    coils.write(output)

    # Read back
    coils2 = CoilsFile.from_file(output)
    assert coils2.n_coils == coils.n_coils
    assert coils2.periods == coils.periods


def test_coilsfile_bounding_box(simple_coils_file):
    """Test bounding box calculation."""
    coils = CoilsFile.from_file(simple_coils_file)
    bbox_min, bbox_max = coils.bounding_box()

    assert bbox_min[0] == 0.0
    assert bbox_max[0] == 0.2
    assert bbox_min[1] == -1.0
    assert bbox_max[1] == 1.0


def test_coilsfile_get_all_points(simple_coils_file):
    """Test get_all_points."""
    coils = CoilsFile.from_file(simple_coils_file)
    pts = coils.get_all_points()

    assert pts.shape == (6, 3)
