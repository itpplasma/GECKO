"""
GECKO - Geometry, Engineering, and CAD for Kinetic modeling and Optimization

A Python library for CAD-related tasks in fusion plasma physics.
"""

__version__ = "0.1.0"

from .coils import (
    CoilsFile,
    Filament,
    coil_to_tube_mesh,
    coils_to_mesh,
    plot_coils_3d,
)
from .collision import (
    CollisionPair,
    CollisionResult,
    coils_to_mesh_with_collisions,
    create_collision_markers,
    detect_coil_collisions,
)

__all__ = [
    "CoilsFile",
    "Filament",
    "coil_to_tube_mesh",
    "coils_to_mesh",
    "plot_coils_3d",
    "CollisionPair",
    "CollisionResult",
    "detect_coil_collisions",
    "create_collision_markers",
    "coils_to_mesh_with_collisions",
]
