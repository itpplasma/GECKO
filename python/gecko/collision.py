"""
Collision detection for fusion device geometry.

This module provides CAD-based collision detection between mesh components,
particularly for detecting overlapping coil winding packs.

Algorithms
----------
The collision detection uses trimesh :cite:p:`trimesh` for mesh intersection
tests. For boolean mesh operations (computing intersection volumes), the Blender
backend is preferred when available. When Blender is not available, a fallback
uses bounding box pre-filtering followed by surface sampling with containment
tests.

For large-scale collision queries, trimesh can use the Flexible Collision
Library (FCL) :cite:p:`pan2012fcl` as a backend for efficient broad-phase
collision detection.

References
----------
.. bibliography::
    :filter: key in {"trimesh", "pan2012fcl"}
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from libneo.coils import CoilsFile

from .coils import coil_to_tube_mesh


@dataclass
class CollisionPair:
    """Result of collision detection between two coils.

    Attributes
    ----------
    coil_i : int
        Index of first coil.
    coil_j : int
        Index of second coil.
    intersection_mesh : trimesh.Trimesh | None
        Mesh of the intersection region, or None if no intersection computed.
    contact_points : np.ndarray
        (N, 3) array of contact/intersection points.
    """

    coil_i: int
    coil_j: int
    intersection_mesh: "trimesh.Trimesh | None"
    contact_points: np.ndarray


@dataclass
class CollisionResult:
    """Result of coil collision analysis.

    Attributes
    ----------
    collisions : list[CollisionPair]
        List of detected collision pairs.
    n_coils : int
        Total number of coils analyzed.
    collision_mesh : trimesh.Trimesh | None
        Combined mesh of all collision regions for visualization.
    """

    collisions: list[CollisionPair]
    n_coils: int
    collision_mesh: "trimesh.Trimesh | None"

    @property
    def has_collisions(self) -> bool:
        """Return True if any collisions were detected."""
        return len(self.collisions) > 0

    @property
    def collision_pairs(self) -> list[tuple[int, int]]:
        """Return list of (i, j) index pairs with collisions."""
        return [(c.coil_i, c.coil_j) for c in self.collisions]


def detect_coil_collisions(
    coils: CoilsFile,
    *,
    width: float = 0.4,
    height: float = 0.8,
) -> CollisionResult:
    """
    Detect collisions between coil winding packs.

    Checks all pairs of coils for geometric intersection of their
    rectangular cross-section tube meshes.

    Parameters
    ----------
    coils : CoilsFile
        CoilsFile object containing coil geometry.
    width : float, default=0.4
        Cross-section width in meters.
    height : float, default=0.8
        Cross-section height in meters.

    Returns
    -------
    CollisionResult
        Result containing collision information and visualization mesh.

    Examples
    --------
    >>> from libneo.coils import CoilsFile
    >>> from gecko.collision import detect_coil_collisions
    >>> coils = CoilsFile.from_file("coils.dat")
    >>> result = detect_coil_collisions(coils, width=0.4, height=0.8)
    >>> if result.has_collisions:
    ...     print(f"Found {len(result.collisions)} collision pairs")
    """
    import trimesh

    # Generate individual coil meshes
    coil_meshes = []
    for fil in coils.filaments:
        if fil.coords.shape[0] < 2:
            coil_meshes.append(None)
            continue
        mesh = coil_to_tube_mesh(fil.coords, width=width, height=height)
        coil_meshes.append(mesh)

    n_coils = len(coil_meshes)
    collisions = []
    intersection_meshes = []

    # Check all pairs for collisions
    for i in range(n_coils):
        if coil_meshes[i] is None:
            continue
        for j in range(i + 1, n_coils):
            if coil_meshes[j] is None:
                continue

            mesh_i = coil_meshes[i]
            mesh_j = coil_meshes[j]

            # Check for intersection using boolean operations
            try:
                intersection = mesh_i.intersection(mesh_j, engine="blender")
                if intersection is not None and len(intersection.vertices) > 0:
                    # Valid intersection found
                    contact_pts = intersection.vertices
                    collisions.append(
                        CollisionPair(
                            coil_i=i,
                            coil_j=j,
                            intersection_mesh=intersection,
                            contact_points=contact_pts,
                        )
                    )
                    intersection_meshes.append(intersection)
            except Exception:
                # Blender not available, fall back to bounding box + sampling
                if mesh_i.bounding_box.intersection(mesh_j.bounding_box) is not None:
                    # Bounding boxes overlap - check more carefully with ray tests
                    contact_pts = _find_contact_points(mesh_i, mesh_j)
                    if len(contact_pts) > 0:
                        collisions.append(
                            CollisionPair(
                                coil_i=i,
                                coil_j=j,
                                intersection_mesh=None,
                                contact_points=contact_pts,
                            )
                        )

    # Combine intersection meshes for visualization
    collision_mesh = None
    if intersection_meshes:
        collision_mesh = trimesh.util.concatenate(intersection_meshes)

    return CollisionResult(
        collisions=collisions,
        n_coils=n_coils,
        collision_mesh=collision_mesh,
    )


def _find_contact_points(
    mesh_a: "trimesh.Trimesh",
    mesh_b: "trimesh.Trimesh",
    n_samples: int = 1000,
) -> np.ndarray:
    """
    Find approximate contact points between two meshes using containment tests.

    Parameters
    ----------
    mesh_a, mesh_b : trimesh.Trimesh
        Meshes to test for intersection.
    n_samples : int
        Number of sample points to test.

    Returns
    -------
    np.ndarray
        (N, 3) array of approximate contact points.
    """
    # Sample points on surface of mesh_a
    points_a, _ = mesh_a.sample(n_samples, return_index=True)

    # Check which points are inside mesh_b
    inside_b = mesh_b.contains(points_a)

    # Sample points on surface of mesh_b
    points_b, _ = mesh_b.sample(n_samples, return_index=True)

    # Check which points are inside mesh_a
    inside_a = mesh_a.contains(points_b)

    # Combine contact points
    contact_from_a = points_a[inside_b]
    contact_from_b = points_b[inside_a]

    if len(contact_from_a) > 0 and len(contact_from_b) > 0:
        return np.vstack([contact_from_a, contact_from_b])
    elif len(contact_from_a) > 0:
        return contact_from_a
    elif len(contact_from_b) > 0:
        return contact_from_b
    else:
        return np.empty((0, 3))


def create_collision_markers(
    result: CollisionResult,
    marker_radius: float = 0.05,
) -> "trimesh.Trimesh":
    """
    Create sphere markers at collision contact points for visualization.

    Use this when boolean intersection is not available (no Blender backend).

    Parameters
    ----------
    result : CollisionResult
        Result from detect_coil_collisions.
    marker_radius : float, default=0.05
        Radius of marker spheres in meters.

    Returns
    -------
    trimesh.Trimesh
        Combined mesh of all marker spheres.
    """
    import trimesh

    markers = []
    for collision in result.collisions:
        for pt in collision.contact_points:
            sphere = trimesh.creation.icosphere(subdivisions=2, radius=marker_radius)
            sphere.apply_translation(pt)
            markers.append(sphere)

    if not markers:
        return trimesh.Trimesh()

    return trimesh.util.concatenate(markers)


def coils_to_mesh_with_collisions(
    coils: CoilsFile,
    *,
    width: float = 0.4,
    height: float = 0.8,
    collision_result: CollisionResult | None = None,
) -> tuple["trimesh.Trimesh", "trimesh.Trimesh | None", set[int]]:
    """
    Generate coil meshes with collision detection, separating colliding coils.

    Parameters
    ----------
    coils : CoilsFile
        CoilsFile object containing coil geometry.
    width : float, default=0.4
        Cross-section width in meters.
    height : float, default=0.8
        Cross-section height in meters.
    collision_result : CollisionResult, optional
        Pre-computed collision result. If None, collision detection is run.

    Returns
    -------
    ok_mesh : trimesh.Trimesh
        Combined mesh of non-colliding coils (render in normal color).
    colliding_mesh : trimesh.Trimesh | None
        Combined mesh of colliding coils (render in warning color), or None.
    colliding_indices : set[int]
        Set of coil indices involved in collisions.
    """
    import trimesh

    from .coils import coil_to_tube_mesh

    if collision_result is None:
        collision_result = detect_coil_collisions(coils, width=width, height=height)

    # Find all coil indices involved in collisions
    colliding_indices = set()
    for c in collision_result.collisions:
        colliding_indices.add(c.coil_i)
        colliding_indices.add(c.coil_j)

    ok_meshes = []
    colliding_meshes = []

    for i, fil in enumerate(coils.filaments):
        if fil.coords.shape[0] < 2:
            continue
        mesh = coil_to_tube_mesh(fil.coords, width=width, height=height)
        if i in colliding_indices:
            colliding_meshes.append(mesh)
        else:
            ok_meshes.append(mesh)

    ok_mesh = trimesh.util.concatenate(ok_meshes) if ok_meshes else trimesh.Trimesh()
    colliding_mesh = (
        trimesh.util.concatenate(colliding_meshes) if colliding_meshes else None
    )

    return ok_mesh, colliding_mesh, colliding_indices


__all__ = [
    "CollisionPair",
    "CollisionResult",
    "detect_coil_collisions",
    "create_collision_markers",
    "coils_to_mesh_with_collisions",
]
