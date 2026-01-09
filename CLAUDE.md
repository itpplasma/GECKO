# CLAUDE.md

This file provides guidance for AI assistants working with GECKO.

## Project Overview

GECKO (Geometry, Engineering, and CAD for Kinetic modeling and Optimization) is a Python library for CAD-related tasks in fusion plasma physics. It extracts and manages geometric data needed for plasma transport simulations.

## Build & Test

```bash
# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest python/tests -v

# Run tests with visual output (generates plots)
pytest python/tests -v --visual
```

## Code Style

- Use type hints for all public functions
- Follow PEP 8 with 88-character line limit
- Prefer composition over inheritance
- Keep modules focused and under 500 lines
- Use dataclasses for data containers
- Lazy imports for heavy dependencies (trimesh, pyvista)

## Module Structure

- `gecko/stl.py` - STL mesh operations (slicing, boundary extraction)
- `gecko/wall.py` - Wall mesh generation from equilibrium data
- `gecko/coils.py` - Coil file parsing and manipulation
- `gecko/viewer.py` - Interactive 3D visualization
- `gecko/cli.py` - Command-line interface

## Dependencies

Core (always available):
- numpy

Optional (lazy-loaded):
- trimesh: STL mesh operations
- shapely: 2D geometry operations
- pyvista: 3D visualization
- matplotlib: 2D plotting
- netCDF4: NetCDF file I/O

## Testing

- Tests should have visual output options for debugging
- Use pytest fixtures for common test data
- Generate and save plots to verify results visually
