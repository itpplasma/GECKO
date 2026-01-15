"""
Wall mesh generation from VMEC equilibria.

This module provides functions for creating 3D wall meshes from VMEC wout files,
with optional port cutouts. Depends on libneo for VMEC geometry reading.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def wall_mesh_from_wout(
    wout_path: Path | str,
    *,
    s_index: int = -1,
    n_theta: int = 128,
    n_zeta: int = 96,
    use_asym: bool = True,
) -> "trimesh.Trimesh":
    """
    Create a 3D triangle mesh of the VMEC LCFS (or other flux surface).

    Args:
        wout_path: Path to VMEC wout NetCDF file
        s_index: Radial surface index (-1 for LCFS)
        n_theta: Number of poloidal points
        n_zeta: Number of toroidal points
        use_asym: Include asymmetric Fourier terms

    Returns:
        trimesh.Trimesh object
    """
    import trimesh

    from libneo.vmec import VMECGeometry

    geom = VMECGeometry.from_file(str(wout_path))
    ns = geom.rmnc.shape[1]

    if s_index < 0:
        s_index = ns + s_index
    if s_index < 0 or s_index >= ns:
        raise ValueError(f"s_index out of range: {s_index} (ns={ns})")

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)

    xyz = np.zeros((n_zeta, n_theta, 3), dtype=float)
    for i, zet in enumerate(zeta):
        r, z, _ = geom.coords(s_index, theta, float(zet), use_asym=use_asym)
        xyz[i, :, 0] = r * np.cos(zet)
        xyz[i, :, 1] = r * np.sin(zet)
        xyz[i, :, 2] = z

    vertices = xyz.reshape((-1, 3))

    faces: list[tuple[int, int, int]] = []
    for i in range(n_zeta):
        inext = (i + 1) % n_zeta
        for j in range(n_theta):
            jnext = (j + 1) % n_theta
            v00 = i * n_theta + j
            v10 = inext * n_theta + j
            v01 = i * n_theta + jnext
            v11 = inext * n_theta + jnext
            faces.append((v00, v10, v11))
            faces.append((v00, v11, v01))

    return trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces), process=False)


def offset_wall_mesh(
    wout_path: Path | str,
    offset: float,
    *,
    n_theta: int = 128,
    n_zeta: int = 96,
    use_asym: bool = True,
) -> "trimesh.Trimesh":
    """
    Create a wall mesh offset from the VMEC LCFS.

    Uses shapely buffer for robust offset on non-convex cross-sections.

    Args:
        wout_path: Path to VMEC wout NetCDF file
        offset: Offset distance (positive = outward) in same units as wout
        n_theta: Number of poloidal points
        n_zeta: Number of toroidal points
        use_asym: Include asymmetric Fourier terms

    Returns:
        trimesh.Trimesh object
    """
    import trimesh
    from shapely.geometry import Polygon

    from libneo.vmec import VMECGeometry

    geom = VMECGeometry.from_file(str(wout_path))
    ns = geom.rmnc.shape[1]
    s_index = ns - 1

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)

    xyz = np.zeros((n_zeta, n_theta, 3), dtype=float)
    for i, zet in enumerate(zeta):
        r_lcfs, z_lcfs, _ = geom.coords(s_index, theta, float(zet), use_asym=use_asym)

        poly = Polygon(np.column_stack([r_lcfs, z_lcfs]))
        offset_poly = poly.buffer(offset, join_style="round")

        if hasattr(offset_poly, "exterior"):
            coords = np.array(offset_poly.exterior.coords)
        else:
            raise RuntimeError("buffer produced invalid geometry")

        r_wall, z_wall = _resample_polygon(coords[:, 0], coords[:, 1], n_theta)

        xyz[i, :, 0] = r_wall * np.cos(zet)
        xyz[i, :, 1] = r_wall * np.sin(zet)
        xyz[i, :, 2] = z_wall

    vertices = xyz.reshape((-1, 3))

    faces: list[tuple[int, int, int]] = []
    for i in range(n_zeta):
        inext = (i + 1) % n_zeta
        for j in range(n_theta):
            jnext = (j + 1) % n_theta
            v00 = i * n_theta + j
            v10 = inext * n_theta + j
            v01 = i * n_theta + jnext
            v11 = inext * n_theta + jnext
            faces.append((v00, v10, v11))
            faces.append((v00, v11, v01))

    return trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces), process=False)


def _resample_polygon(
    r: np.ndarray, z: np.ndarray, n: int
) -> tuple[np.ndarray, np.ndarray]:
    """Resample a closed polygon to n evenly-spaced points."""
    pts = np.column_stack([r, z])
    if np.linalg.norm(pts[0] - pts[-1]) < 1e-10:
        pts = pts[:-1]

    ds = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    ds = np.append(ds, np.linalg.norm(pts[0] - pts[-1]))
    s = np.concatenate(([0.0], np.cumsum(ds)))
    total = s[-1]

    s_new = np.linspace(0.0, total, n, endpoint=False)
    r_new = np.interp(s_new, s, np.r_[pts[:, 0], pts[0, 0]])
    z_new = np.interp(s_new, s, np.r_[pts[:, 1], pts[0, 1]])

    return r_new, z_new


@dataclass(frozen=True)
class PortSpec:
    """Specification for a cylindrical port cutout."""

    phi: float
    z: float
    radius: float
    length: float
    inboard: bool = False


def _wrap_angle_pi(angle: np.ndarray) -> np.ndarray:
    return (angle + np.pi) % (2.0 * np.pi) - np.pi


def _find_surface_point(
    mesh: "trimesh.Trimesh",
    phi_target: float,
    z_target: float,
    *,
    phi_tol: float = 0.15,
    z_tol: float = 0.3,
    inboard: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Find a point on the mesh surface and outward direction."""
    centers = mesh.triangles_center
    phi = np.arctan2(centers[:, 1], centers[:, 0])
    z = centers[:, 2]
    r = np.sqrt(centers[:, 0] ** 2 + centers[:, 1] ** 2)

    dphi = np.abs(_wrap_angle_pi(phi - float(phi_target)))
    dz = np.abs(z - float(z_target))
    mask = (dphi <= phi_tol) & (dz <= z_tol)
    if not np.any(mask):
        raise RuntimeError("no candidate faces for requested (phi,z)")

    if inboard:
        idx = int(np.argmin(np.where(mask, r, np.inf)))
        sign = -1.0
    else:
        idx = int(np.argmax(np.where(mask, r, -np.inf)))
        sign = 1.0

    center = centers[idx].copy()
    direction = sign * np.array(
        [np.cos(phi_target), np.sin(phi_target), 0.0], dtype=float
    )
    nrm = float(np.linalg.norm(direction))
    if nrm == 0.0:
        raise RuntimeError("degenerate port direction")
    direction /= nrm
    return center, direction


def cut_cylindrical_ports(
    mesh: "trimesh.Trimesh", ports: list[PortSpec]
) -> "trimesh.Trimesh":
    """
    Cut cylindrical ports from a wall mesh.

    Args:
        mesh: Input trimesh object
        ports: List of PortSpec defining port locations

    Returns:
        New trimesh with port faces removed
    """
    import trimesh

    centers = mesh.triangles_center
    keep = np.ones((mesh.faces.shape[0],), dtype=bool)

    for p in ports:
        center, direction = _find_surface_point(mesh, p.phi, p.z, inboard=p.inboard)
        v = centers - center[None, :]
        t = v @ direction
        radial = np.linalg.norm(v - t[:, None] * direction[None, :], axis=1)
        half_len = 0.5 * float(p.length)
        in_cyl = (np.abs(t) <= half_len) & (radial <= float(p.radius))
        keep &= ~in_cyl

    faces = mesh.faces[keep]
    return trimesh.Trimesh(vertices=mesh.vertices.copy(), faces=faces, process=False)


def wall_mesh_from_chartmap(
    chartmap_path: Path | str,
    *,
    rho_index: int = -1,
    scale_to_m: float = 0.01,
) -> "trimesh.Trimesh":
    """
    Create a 3D triangle mesh from a chartmap boundary surface.

    This creates an STL wall that exactly matches the chartmap coordinate
    boundary at the specified rho index. Use this instead of offset_wall_mesh
    when the chartmap was generated with map2disc or vmec_extended, as the
    boundary shape may differ from a simple VMEC LCFS offset.

    Args:
        chartmap_path: Path to chartmap NetCDF file
        rho_index: Radial index for boundary (-1 for outermost = rho=1)
        scale_to_m: Scale factor to convert chartmap units to meters
                    (0.01 for cm, 1.0 for m)

    Returns:
        trimesh.Trimesh object with vertices in meters
    """
    import trimesh
    from netCDF4 import Dataset

    with Dataset(str(chartmap_path), "r") as ds:
        theta = ds.variables["theta"][:]
        zeta = ds.variables["zeta"][:]
        nfp = int(ds.variables["num_field_periods"][...])

        x = ds.variables["x"][:, :, rho_index]  # (nzeta, ntheta)
        y = ds.variables["y"][:, :, rho_index]
        z = ds.variables["z"][:, :, rho_index]

    nzeta, ntheta = x.shape

    # Replicate for all field periods
    all_x, all_y, all_z = [], [], []
    for fp in range(nfp):
        angle = 2.0 * np.pi * fp / nfp
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        x_rot = x * cos_a - y * sin_a
        y_rot = x * sin_a + y * cos_a
        all_x.append(x_rot)
        all_y.append(y_rot)
        all_z.append(z)

    x_full = np.concatenate(all_x, axis=0)
    y_full = np.concatenate(all_y, axis=0)
    z_full = np.concatenate(all_z, axis=0)
    nzeta_total = nzeta * nfp

    # Assume last theta wraps to first (theta=2pi == theta=0)
    ntheta_unique = ntheta - 1 if ntheta > 1 else ntheta

    # Build vertices
    vertices = np.zeros((nzeta_total * ntheta_unique, 3))
    for iz in range(nzeta_total):
        for it in range(ntheta_unique):
            idx = iz * ntheta_unique + it
            vertices[idx] = [x_full[iz, it], y_full[iz, it], z_full[iz, it]]

    # Convert to meters
    vertices *= scale_to_m

    # Build faces (wrapping both theta and zeta)
    faces = []
    for iz in range(nzeta_total):
        iz_next = (iz + 1) % nzeta_total
        for it in range(ntheta_unique):
            it_next = (it + 1) % ntheta_unique

            i00 = iz * ntheta_unique + it
            i10 = iz_next * ntheta_unique + it
            i01 = iz * ntheta_unique + it_next
            i11 = iz_next * ntheta_unique + it_next

            faces.append([i00, i10, i01])
            faces.append([i10, i11, i01])

    mesh = trimesh.Trimesh(
        vertices=vertices, faces=np.asarray(faces), process=False
    )
    mesh.fix_normals()

    return mesh
