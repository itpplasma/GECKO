"""
Wall mesh generation from VMEC equilibria.

This module provides functions for creating 3D wall meshes from VMEC wout files,
with optional port cutouts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _read_vmec_boundary(wout_path: Path | str) -> dict:
    """
    Read VMEC boundary Fourier coefficients from a wout NetCDF file.

    Returns:
        Dictionary with keys: rmnc, zmns, rmns, zmnc, xm, xn, nfp, ns
    """
    import netCDF4 as nc

    wout_path = Path(wout_path)
    with nc.Dataset(wout_path, "r") as ds:
        rmnc = np.asarray(ds.variables["rmnc"][:], dtype=float)
        zmns = np.asarray(ds.variables["zmns"][:], dtype=float)
        xm = np.asarray(ds.variables["xm"][:], dtype=float)
        xn = np.asarray(ds.variables["xn"][:], dtype=float)
        nfp = int(ds.variables["nfp"][...])
        ns = int(ds.variables["ns"][...])

        rmns = None
        zmnc = None
        if "rmns" in ds.variables:
            rmns = np.asarray(ds.variables["rmns"][:], dtype=float)
        if "zmnc" in ds.variables:
            zmnc = np.asarray(ds.variables["zmnc"][:], dtype=float)

    return {
        "rmnc": rmnc,
        "zmns": zmns,
        "rmns": rmns,
        "zmnc": zmnc,
        "xm": xm,
        "xn": xn,
        "nfp": nfp,
        "ns": ns,
    }


def _vmec_coords(
    vmec_data: dict,
    s_index: int,
    theta: np.ndarray,
    zeta: float,
    *,
    use_asym: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute R, Z coordinates from VMEC Fourier coefficients.

    Args:
        vmec_data: Dictionary from _read_vmec_boundary
        s_index: Radial surface index
        theta: Poloidal angles (1D array)
        zeta: Single toroidal angle
        use_asym: Include asymmetric terms if available

    Returns:
        R, Z arrays of same shape as theta
    """
    rmnc = vmec_data["rmnc"]
    zmns = vmec_data["zmns"]
    xm = vmec_data["xm"]
    xn = vmec_data["xn"]
    rmns = vmec_data["rmns"]
    zmnc = vmec_data["zmnc"]

    theta = np.asarray(theta, dtype=float)
    angle = xm[:, None] * theta[None, :] - xn[:, None] * zeta

    r = np.sum(rmnc[:, s_index, None] * np.cos(angle), axis=0)
    z = np.sum(zmns[:, s_index, None] * np.sin(angle), axis=0)

    if use_asym and rmns is not None and zmnc is not None:
        r += np.sum(rmns[:, s_index, None] * np.sin(angle), axis=0)
        z += np.sum(zmnc[:, s_index, None] * np.cos(angle), axis=0)

    return r, z


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

    vmec_data = _read_vmec_boundary(wout_path)
    ns = vmec_data["ns"]

    if s_index < 0:
        s_index = ns + s_index
    if s_index < 0 or s_index >= ns:
        raise ValueError(f"s_index out of range: {s_index} (ns={ns})")

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)

    xyz = np.zeros((n_zeta, n_theta, 3), dtype=float)
    for i, zet in enumerate(zeta):
        r, z = _vmec_coords(vmec_data, s_index, theta, float(zet), use_asym=use_asym)
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

    vmec_data = _read_vmec_boundary(wout_path)
    ns = vmec_data["ns"]
    s_index = ns - 1

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)

    xyz = np.zeros((n_zeta, n_theta, 3), dtype=float)
    for i, zet in enumerate(zeta):
        r_lcfs, z_lcfs = _vmec_coords(
            vmec_data, s_index, theta, float(zet), use_asym=use_asym
        )

        poly = Polygon(np.column_stack([r_lcfs, z_lcfs]))
        offset_poly = poly.buffer(offset, join_style=1)

        if hasattr(offset_poly, "exterior"):
            coords = np.array(offset_poly.exterior.coords)
        else:
            raise RuntimeError("buffer produced invalid geometry")

        # Resample to n_theta points
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


def _resample_polygon(r: np.ndarray, z: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
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
