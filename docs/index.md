# GECKO

**G**eometry, **E**ngineering, and **C**AD for **K**inetic modeling and **O**ptimization

GECKO provides CAD and visualization tools for fusion device geometry, building on
[libneo](https://github.com/itpplasma/libneo) for core physics functionality.

## Features

- **Wall mesh generation**: Create STL meshes from VMEC equilibria with offset support
- **Coil visualization**: 3D tube meshes with rectangular cross-sections
- **Interactive viewer**: PyVista-based 3D visualization

## Installation

```bash
pip install gecko
```

For development:

```bash
git clone git@github.com:itpplasma/GECKO.git
cd GECKO
pip install -e ".[dev]"
```

## Quick Start

```python
from gecko.viewer import quick_view

# Visualize coils and wall
quick_view("wall.stl", "coils.scaled")
```

## Documentation

```{toctree}
:maxdepth: 2

methods
api/index
references
```

## References

```{bibliography}
:all:
```
