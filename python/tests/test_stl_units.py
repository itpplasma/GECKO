from __future__ import annotations

from pathlib import Path


def test_convert_stl_units_mm_to_m(tmp_path: Path) -> None:
    import numpy as np
    import trimesh

    from gecko.stl import convert_stl_units

    mesh_m = trimesh.creation.box(extents=(2.0, 4.0, 6.0))
    mesh_mm = mesh_m.copy()
    mesh_mm.apply_scale(1000.0)

    inp = tmp_path / "box_mm.stl"
    out = tmp_path / "box_m.stl"
    mesh_mm.export(str(inp))

    convert_stl_units(inp, out, from_units="mm", to_units="m")

    loaded = trimesh.load_mesh(out, process=False)
    assert isinstance(loaded, trimesh.Trimesh)

    extents = loaded.bounding_box.extents
    assert np.allclose(extents, [2.0, 4.0, 6.0], rtol=0.0, atol=1.0e-12)

