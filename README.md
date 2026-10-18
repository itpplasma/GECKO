# GECKO

**Geometry, Engineering, and CAD for Kinetic modeling and Optimization**

A Python library for handling CAD-related tasks in fusion plasma physics:
- STL mesh slicing and boundary extraction
- Wall geometry generation from VMEC equilibria
- Coil file parsing and visualization
- Interactive 3D visualization with PyVista

## Installation

```bash
pip install -e .

# With STL support (trimesh, shapely)
pip install -e ".[stl]"

# With visualization (pyvista, matplotlib)
pip install -e ".[vis]"

# Full development install
pip install -e ".[dev]"
```

## Features

### Wall Generation from VMEC

Generate STL wall meshes from VMEC equilibrium files:

```python
from gecko.wall import wall_mesh_from_wout

mesh = wall_mesh_from_wout("wout.nc", n_theta=128, n_zeta=96)
mesh.export("wall.stl")
```

### STL Boundary Extraction

Extract RZ cross-sections from STL wall files:

```python
from gecko.stl import extract_boundary_slices

slices = extract_boundary_slices("wall.stl", n_phi=32)
```

### Coil File Handling

Parse and visualize STELLOPT/MAKEGRID coil files:

```python
from gecko.coils import CoilsFile

coils = CoilsFile.from_file("coils.c09r00")
fig, ax = coils.plot(show=True)
```

### Interactive 3D Visualization

```python
from gecko.viewer import Viewer

viewer = Viewer()
viewer.add_wall("wall.stl")
viewer.add_coils("coils.scaled", cross_section=(0.4, 0.8))
viewer.show()
```

### STL Units (mm vs m)

STL files do not store units. GECKO assumes **meters** by default.

For engineer CAD STLs (often **mm**), either:

```bash
# 1) Convert once (recommended)
python3 -m gecko.cli stl-convert wall_mm.stl wall_m.stl --from-units mm --to-units m

# 2) Or scale for visualization only
python3 -m gecko.cli view --stl-units mm wall_mm.stl
```

## License

MIT License

## Public tests and STL provenance

CI builds [libneo](https://github.com/itpplasma/libneo) from a pinned public
commit, installs `.[dev]`, and runs `python -m pytest python/tests -v -ra`.
The system and Python build dependencies are listed in
[the workflow](.github/workflows/tests.yml). libneo is currently installed from
source rather than the Python package index. The build uses at most eight jobs.
For example, after installing the workflow's system dependencies:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy netCDF4 scikit-build-core f90wrap==0.3.0
git clone https://github.com/itpplasma/libneo.git .ci-deps/libneo
git -C .ci-deps/libneo checkout f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
CMAKE_BUILD_PARALLEL_LEVEL=8 FC=gfortran python -m pip install \
  --no-build-isolation .ci-deps/libneo -Ccmake.define.LIBNEO_BUILD_TESTING=OFF
python -m pip install -e '.[dev]'
python -m pytest python/tests -v -ra
```

The vendored [ALPES fixture](python/tests/data/alpes/README.md) is attributed to
The ALPES Team and licensed CC BY 4.0, separately from GECKO's MIT code.
Tests check closed, outward-oriented offset walls, RZ clearance from a densely
sampled LCFS, coil tube profiles, and separated/overlapping coil pairs.
Private geometry tests skip by default. To run them, provide
`GECKO_PRIVATE_WOUT` and `GECKO_PRIVATE_COILS` and pass `--private-data`.
The private viewer test also requires a working off-screen rendering backend.

`gecko wall-from-vmec`, `gecko stl-convert`, and
`gecko.stl.export_mesh_to_stl(mesh, path)` write a `<filename>.stl.json` sidecar
with the GECKO package version, `git describe --tags --always --dirty`, and the
STL SHA-256 checksum. Keep the sidecar with meshes used in publications. Wheels
retain their build revision; source archives built without Git record a null
revision. Direct third-party calls such as `mesh.export(...)` bypass this writer.
