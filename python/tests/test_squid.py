"""
Test GECKO with SQUID 5C configuration.

This test requires SQUID data files at:
- ~/data/SQUID/VMEC/Proxima2024/5C_6NKvZumZYmPEDi3rmLRd/VmecFreeBoundary/wout_free.nc
- ~/data/SQUID/VMEC/Proxima2024/5C_6NKvZumZYmPEDi3rmLRd/ScaleCoils/coils.scaled

Run with: pytest python/tests/test_squid.py -v --visual

To generate figures:
    python python/tests/test_squid.py
"""

from pathlib import Path

import pytest

SQUID_WOUT = Path.home() / "data/SQUID/VMEC/Proxima2024/5C_6NKvZumZYmPEDi3rmLRd/VmecFreeBoundary/wout_free.nc"
SQUID_COILS = Path.home() / "data/SQUID/VMEC/Proxima2024/5C_6NKvZumZYmPEDi3rmLRd/ScaleCoils/coils.scaled"
OUTPUT_DIR = Path("/tmp/gecko_squid")


@pytest.fixture
def squid_data():
    """Skip test if SQUID data not available."""
    if not SQUID_WOUT.exists() or not SQUID_COILS.exists():
        pytest.skip("SQUID data not available")
    return SQUID_WOUT, SQUID_COILS


def test_squid_wall_mesh(squid_data):
    """Test creating wall mesh from SQUID VMEC."""
    pytest.importorskip("trimesh")

    from gecko.wall import wall_mesh_from_wout

    wout, _ = squid_data
    mesh = wall_mesh_from_wout(wout, n_theta=64, n_zeta=48)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_squid_offset_wall(squid_data):
    """Test creating offset wall (25cm from LCFS)."""
    pytest.importorskip("trimesh")
    pytest.importorskip("shapely")

    from gecko.wall import offset_wall_mesh

    wout, _ = squid_data
    mesh = offset_wall_mesh(wout, offset=0.25, n_theta=64, n_zeta=48)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mesh.export(str(OUTPUT_DIR / "squid_wall_25cm.stl"))


def test_squid_coils_mesh(squid_data):
    """Test creating coil tube meshes."""
    pytest.importorskip("trimesh")

    from libneo.coils import CoilsFile

    from gecko.coils import coils_to_mesh

    _, coils_path = squid_data
    coils = CoilsFile.from_file(str(coils_path))
    mesh = coils_to_mesh(coils, width=0.4, height=0.8)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mesh.export(str(OUTPUT_DIR / "squid_coils.stl"))


def test_squid_viewer(squid_data):
    """Test full SQUID visualization (requires display)."""
    pytest.importorskip("pyvista")

    from gecko.viewer import Viewer

    wout, coils_path = squid_data

    viewer = Viewer()
    viewer.add_vmec_lcfs(wout, color="red", opacity=0.3, n_theta=64, n_zeta=48)
    viewer.add_coils(coils_path, cross_section=(0.4, 0.8), color="copper")
    viewer.add_axes()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    viewer.screenshot(OUTPUT_DIR / "squid_visualization.png", window_size=(1920, 1080))
    viewer.close()


def main():
    """Generate SQUID visualization figures."""
    import pyvista as pv
    import trimesh

    from libneo.coils import CoilsFile

    from gecko.coils import coils_to_mesh
    from gecko.viewer import Viewer
    from gecko.wall import offset_wall_mesh, wall_mesh_from_wout

    if not SQUID_WOUT.exists() or not SQUID_COILS.exists():
        print("SQUID data not found. Expected:")
        print(f"  {SQUID_WOUT}")
        print(f"  {SQUID_COILS}")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Creating LCFS mesh...")
    lcfs_mesh = wall_mesh_from_wout(SQUID_WOUT, n_theta=128, n_zeta=96)
    lcfs_mesh.export(str(OUTPUT_DIR / "squid_lcfs.stl"))
    print(f"  Saved: {OUTPUT_DIR / 'squid_lcfs.stl'}")

    print("Creating wall mesh (25cm offset)...")
    wall_mesh = offset_wall_mesh(SQUID_WOUT, offset=0.25, n_theta=128, n_zeta=96)
    wall_mesh.export(str(OUTPUT_DIR / "squid_wall_25cm.stl"))
    print(f"  Saved: {OUTPUT_DIR / 'squid_wall_25cm.stl'}")

    print("Creating coil meshes (SPARC-like: 40cm x 80cm)...")
    coils = CoilsFile.from_file(str(SQUID_COILS))
    coil_mesh = coils_to_mesh(coils, width=0.4, height=0.8)
    coil_mesh.export(str(OUTPUT_DIR / "squid_coils.stl"))
    print(f"  Saved: {OUTPUT_DIR / 'squid_coils.stl'}")
    print(f"  Number of coils: {len(coils.filaments)}")
    print(f"  Field periods: {coils.periods}")

    print("Creating visualization (off-screen rendering)...")
    viewer = Viewer(off_screen=True)
    viewer.add_wall(lcfs_mesh, color="red", opacity=0.4)
    viewer.add_wall(wall_mesh, color="lightgray", opacity=0.3)
    viewer.add_coils(coils, cross_section=(0.4, 0.8), color="#B87333")
    viewer.add_axes()

    # Save screenshot (off-screen mode)
    viewer.screenshot(OUTPUT_DIR / "squid_full.png", window_size=(1920, 1080))
    print(f"  Saved: {OUTPUT_DIR / 'squid_full.png'}")

    # Close viewer (no display available in headless mode)
    viewer.close()

    print("\nTo view interactively, run:")
    print(f"  gecko view {OUTPUT_DIR / 'squid_lcfs.stl'} {OUTPUT_DIR / 'squid_wall_25cm.stl'} {OUTPUT_DIR / 'squid_coils.stl'}")
    print("\nOr from Python:")
    print("  from gecko.viewer import quick_view")
    print(f"  quick_view('{OUTPUT_DIR / 'squid_lcfs.stl'}', '{OUTPUT_DIR / 'squid_coils.stl'}')")


if __name__ == "__main__":
    main()
