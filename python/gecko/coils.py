"""
Coil visualization and mesh generation for fusion device geometry.

This module extends libneo.coils with 3D visualization capabilities:
- Rectangular cross-section tube meshes for coils
- PyVista/trimesh integration

Core coil parsing is provided by libneo.coils.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Re-export from libneo for convenience
from libneo.coils import CoilsFile, Filament


def coil_to_tube_mesh(
    coords: np.ndarray,
    *,
    width: float,
    height: float,
    close_threshold: float = 0.1,
) -> "trimesh.Trimesh":
    """
    Create a 3D tube mesh from coil centerline with rectangular cross-section.

    Args:
        coords: (N, 3) coil centerline points in meters
        width: Cross-section width (horizontal) in meters
        height: Cross-section height (vertical) in meters
        close_threshold: If first-last point distance < threshold, close the loop

    Returns:
        trimesh.Trimesh object
    """
    import trimesh

    if coords.shape[0] < 2:
        raise ValueError("Need at least 2 points for tube mesh")

    # Check if coil is closed (first point close to last point)
    gap = np.linalg.norm(coords[0] - coords[-1])
    is_closed = gap < close_threshold

    # For closed coils, exclude the duplicate last point
    if is_closed:
        coords = coords[:-1]

    n_pts = coords.shape[0]
    if n_pts < 2:
        raise ValueError("Need at least 2 points for tube mesh after removing duplicate")

    # Compute tangent vectors
    tangent = np.zeros_like(coords)
    if is_closed:
        # For closed coils, use periodic tangents (wrap around)
        tangent[0] = coords[1] - coords[-1]
        tangent[-1] = coords[0] - coords[-2]
    else:
        tangent[0] = coords[1] - coords[0]
        tangent[-1] = coords[-1] - coords[-2]
    tangent[1:-1] = coords[2:] - coords[:-2]
    norms = np.linalg.norm(tangent, axis=1, keepdims=True)
    norms[norms == 0] = 1e-12
    tangent = tangent / norms

    # Create cross-section profile (rectangular)
    half_w = width / 2
    half_h = height / 2
    profile = np.array([
        [-half_w, -half_h],
        [half_w, -half_h],
        [half_w, half_h],
        [-half_w, half_h],
    ])
    n_profile = len(profile)

    vertices_list = []
    up = np.array([0.0, 0.0, 1.0])

    for center, t in zip(coords, tangent):
        # Compute local frame using cross product
        if abs(np.dot(t, up)) > 0.99:
            side = np.cross(t, np.array([1.0, 0.0, 0.0]))
        else:
            side = np.cross(t, up)
        side_norm = np.linalg.norm(side)
        if side_norm > 1e-12:
            side = side / side_norm
        else:
            side = np.array([1.0, 0.0, 0.0])

        up_local = np.cross(side, t)
        up_norm = np.linalg.norm(up_local)
        if up_norm > 1e-12:
            up_local = up_local / up_norm
        else:
            up_local = np.array([0.0, 0.0, 1.0])

        # Transform profile to world coordinates
        for px, py in profile:
            pt = center + px * side + py * up_local
            vertices_list.append(pt)

    vertices = np.array(vertices_list)
    faces = []

    # Create faces connecting adjacent profiles
    for i in range(n_pts - 1):
        base = i * n_profile
        next_base = (i + 1) * n_profile
        for j in range(n_profile):
            j_next = (j + 1) % n_profile
            v0 = base + j
            v1 = base + j_next
            v2 = next_base + j_next
            v3 = next_base + j
            faces.append([v0, v1, v2])
            faces.append([v0, v2, v3])

    if is_closed:
        # Connect last profile back to first profile
        base = (n_pts - 1) * n_profile
        next_base = 0
        for j in range(n_profile):
            j_next = (j + 1) % n_profile
            v0 = base + j
            v1 = base + j_next
            v2 = next_base + j_next
            v3 = next_base + j
            faces.append([v0, v1, v2])
            faces.append([v0, v2, v3])
    else:
        # Cap the open ends
        def add_cap(base_idx: int, reverse: bool = False):
            nonlocal vertices
            cap_verts = vertices[base_idx : base_idx + n_profile]
            center = cap_verts.mean(axis=0)
            center_idx = len(vertices)
            vertices = np.vstack([vertices, center])
            for j in range(n_profile):
                j_next = (j + 1) % n_profile
                if reverse:
                    faces.append([center_idx, base_idx + j_next, base_idx + j])
                else:
                    faces.append([center_idx, base_idx + j, base_idx + j_next])

        add_cap(0, reverse=True)
        add_cap((n_pts - 1) * n_profile, reverse=False)

    return trimesh.Trimesh(vertices=vertices, faces=np.array(faces), process=False)


def coils_to_mesh(
    coils: CoilsFile,
    *,
    width: float = 0.4,
    height: float = 0.8,
) -> "trimesh.Trimesh":
    """
    Convert all coils to a combined 3D mesh with rectangular cross-sections.

    Args:
        coils: CoilsFile object (from libneo.coils or gecko.coils)
        width: Cross-section width in meters
        height: Cross-section height in meters

    Returns:
        Combined trimesh.Trimesh object
    """
    import trimesh

    meshes = []
    for fil in coils.filaments:
        if fil.coords.shape[0] < 2:
            continue
        mesh = coil_to_tube_mesh(fil.coords, width=width, height=height)
        meshes.append(mesh)

    if not meshes:
        return trimesh.Trimesh()

    return trimesh.util.concatenate(meshes)


def plot_coils_3d(
    coils: CoilsFile,
    *,
    cross_section: tuple[float, float] | None = None,
    color: str = "copper",
    show: bool = True,
):
    """
    Plot coils in 3D using PyVista.

    Args:
        coils: CoilsFile object
        cross_section: (width, height) for tube mesh, or None for centerline only
        color: Mesh/line color
        show: Whether to display immediately

    Returns:
        PyVista plotter
    """
    import pyvista as pv

    plotter = pv.Plotter()
    plotter.set_background("white")

    if cross_section is not None:
        width, height = cross_section
        mesh = coils_to_mesh(coils, width=width, height=height)
        pv_mesh = pv.wrap(mesh)
        plotter.add_mesh(pv_mesh, color=color)
    else:
        for fil in coils.filaments:
            if fil.coords.shape[0] < 2:
                continue
            line = pv.lines_from_points(fil.coords)
            plotter.add_mesh(line, color=color, line_width=2.0)

    plotter.add_axes()

    if show:
        plotter.show()

    return plotter


__all__ = [
    "CoilsFile",
    "Filament",
    "coil_to_tube_mesh",
    "coils_to_mesh",
    "plot_coils_3d",
]
