# GECKO

Geometry, Engineering, and CAD for Kinetic modeling and Optimization: a Python
library for CAD-related tasks in fusion plasma physics, built on
[libneo](https://github.com/itpplasma/libneo).

- STL mesh slicing and boundary extraction
- Wall geometry generation from VMEC equilibria
- Coil file parsing and visualization
- Coil filament extraction from STEP winding packs (CadQuery)
- Interactive 3D visualization with PyVista

Status: alpha; interfaces may change.

## Installation

```bash
pip install -e .

# With STL support (trimesh, shapely)
pip install -e ".[stl]"

# With visualization (pyvista, matplotlib)
pip install -e ".[vis]"

# With STEP support (cadquery)
pip install -e ".[step]"

# Full development install
pip install -e ".[dev]"
```

CI installs the `libneo` dependency from source at a pinned commit; see
[docs/testing.md](docs/testing.md).

## Features

### Wall generation from VMEC

Generate STL wall meshes from VMEC equilibrium files:

```python
from gecko.wall import wall_mesh_from_wout

mesh = wall_mesh_from_wout("wout.nc", n_theta=128, n_zeta=96)
mesh.export("wall.stl")
```

### STL boundary extraction

Extract RZ cross-sections from STL wall files:

```python
from gecko.stl import extract_boundary_slices

slices = extract_boundary_slices("wall.stl", n_phi=32)
```

### Coil files

Parse and visualize STELLOPT/MAKEGRID coil files:

```python
from gecko.coils import CoilsFile

coils = CoilsFile.from_file("coils.c09r00")
fig, ax = coils.plot(show=True)
```

### Interactive 3D visualization

```python
from gecko.viewer import Viewer

viewer = Viewer()
viewer.add_wall("wall.stl")
viewer.add_coils("coils.scaled", cross_section=(0.4, 0.8))
viewer.show()
```

### STL units

STL files do not store units. GECKO assumes metres by default. For CAD STL
files in millimetres, either:

```bash
# 1) Convert once (recommended)
python3 -m gecko.cli stl-convert wall_mm.stl wall_m.stl --from-units mm --to-units m

# 2) Or scale for visualization only
python3 -m gecko.cli view --stl-units mm wall_mm.stl
```

## Command line

The `gecko` command (also `python3 -m gecko.cli`) provides `view`,
`wall-from-vmec`, `coils-from-step`, `coils-info`, `coils-plot`,
`stl-convert` and `view-results` (SIMPLE wall hits on an STL wall); see
`gecko --help`.

## Testing

    pip install -e ".[dev,step]"
    python -m pytest python/tests -v -ra

See [docs/testing.md](docs/testing.md) for the full setup used in CI,
including libneo, and for the optional private-data tests.

## Documentation

The Sphinx documentation in `docs/` covers the API, the
[methods](docs/methods.md) used for wall and coil meshes,
[testing](docs/testing.md), [STL provenance](docs/provenance.md) and
[contributing](docs/contributing.md).

## License

MIT, see [LICENSE](LICENSE). The ALPES test fixture in `python/tests/data/alpes/`
is licensed CC BY 4.0 (see its [README](python/tests/data/alpes/README.md)).
