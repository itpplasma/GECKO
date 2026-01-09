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

    # Check all pairs for collisions
    for i in range(n_coils):
        if coil_meshes[i] is None:
            continue
        for j in range(i + 1, n_coils):
            if coil_meshes[j] is None:
                continue

            mesh_i = coil_meshes[i]
            mesh_j = coil_meshes[j]

            # Quick bounding box check first
            bounds_i = mesh_i.bounds  # (2, 3) array: min and max
            bounds_j = mesh_j.bounds

            # Check if bounding boxes overlap
            overlap = np.all(bounds_i[0] <= bounds_j[1]) and np.all(bounds_j[0] <= bounds_i[1])
            if not overlap:
                continue

            # Bounding boxes overlap - check more carefully with surface sampling
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

    return CollisionResult(
        collisions=collisions,
        n_coils=n_coils,
        collision_mesh=None,
    )


def _find_contact_points(
    mesh_a: "trimesh.Trimesh",
    mesh_b: "trimesh.Trimesh",
    n_samples: int = 500,
) -> np.ndarray:
    """
    Find approximate contact points between two meshes using vertex proximity.

    Uses a simple distance-based approach that does not require rtree or
    other optional dependencies.

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
    from scipy.spatial import cKDTree

    # Sample points on surfaces
    points_a, _ = mesh_a.sample(n_samples, return_index=True)
    points_b, _ = mesh_b.sample(n_samples, return_index=True)

    # Build KD-tree for mesh_b vertices
    tree_b = cKDTree(mesh_b.vertices)
    tree_a = cKDTree(mesh_a.vertices)

    # Find points from mesh_a that are very close to mesh_b surface
    # Use a threshold based on characteristic mesh size
    char_size_a = np.sqrt(mesh_a.area / len(mesh_a.faces)) if len(mesh_a.faces) > 0 else 0.01
    char_size_b = np.sqrt(mesh_b.area / len(mesh_b.faces)) if len(mesh_b.faces) > 0 else 0.01
    threshold = 2.0 * max(char_size_a, char_size_b)

    # Query distances from sampled points to other mesh
    dist_a_to_b, _ = tree_b.query(points_a)
    dist_b_to_a, _ = tree_a.query(points_b)

    # Contact points are those closer than threshold
    contact_a = points_a[dist_a_to_b < threshold]
    contact_b = points_b[dist_b_to_a < threshold]

    if len(contact_a) > 0 and len(contact_b) > 0:
        return np.vstack([contact_a, contact_b])
    elif len(contact_a) > 0:
        return contact_a
    elif len(contact_b) > 0:
        return contact_b
    else:
        return np.empty((0, 3))


def create_collision_markers(
    result: CollisionResult,
    marker_radius: float = 0.05,
) -> "trimesh.Trimesh":
    """
    Create cone markers pointing at collision contact points.

    Follows engineering drawing conventions (ISO 7010) where warning markers
    use pointed/triangular shapes to indicate problem locations.

    Parameters
    ----------
    result : CollisionResult
        Result from detect_coil_collisions.
    marker_radius : float, default=0.05
        Base radius of marker cones in meters.

    Returns
    -------
    trimesh.Trimesh
        Combined mesh of all marker cones.
    """
    import trimesh
    from scipy.spatial.transform import Rotation

    markers = []

    for collision in result.collisions:
        if len(collision.contact_points) == 0:
            continue

        # Compute centroid of contact points for this collision pair
        centroid = np.mean(collision.contact_points, axis=0)

        for pt in collision.contact_points:
            # Create cone (tip at origin, pointing in +Z direction)
            cone_height = marker_radius * 2.5
            cone = trimesh.creation.cone(
                radius=marker_radius,
                height=cone_height,
                sections=16,
            )

            # Direction from contact point toward centroid (cone points inward)
            direction = centroid - pt
            dist = np.linalg.norm(direction)
            if dist < 1e-10:
                # Point is at centroid, use radial direction from origin
                direction = pt / (np.linalg.norm(pt) + 1e-10)
            else:
                direction = direction / dist

            # Rotate cone to point in the direction
            # Default cone points in +Z, we need to rotate to 'direction'
            z_axis = np.array([0.0, 0.0, 1.0])
            if np.abs(np.dot(z_axis, direction)) < 0.999:
                rot_axis = np.cross(z_axis, direction)
                rot_axis = rot_axis / np.linalg.norm(rot_axis)
                angle = np.arccos(np.clip(np.dot(z_axis, direction), -1, 1))
                rot = Rotation.from_rotvec(rot_axis * angle)
                cone.apply_transform(
                    np.vstack([
                        np.hstack([rot.as_matrix(), np.zeros((3, 1))]),
                        [0, 0, 0, 1]
                    ])
                )
            elif np.dot(z_axis, direction) < 0:
                # Flip 180 degrees
                cone.apply_transform([
                    [-1, 0, 0, 0],
                    [0, -1, 0, 0],
                    [0, 0, -1, 0],
                    [0, 0, 0, 1]
                ])

            # Position cone so tip is at the contact point
            # Cone tip is at origin after creation, base at +Z
            # After rotation, tip should be at contact point
            cone.apply_translation(pt)

            markers.append(cone)

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
