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

## License

MIT License
