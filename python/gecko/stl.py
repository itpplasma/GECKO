"""
STL mesh operations for fusion geometry.

This module provides functions for:
- Slicing STL meshes at constant toroidal angle
- Extracting and filling RZ boundary contours
- Stitching open contours with C1 bridges
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np


def _unit(vec: np.ndarray) -> np.ndarray:
    nrm = float(np.linalg.norm(vec))
    if nrm == 0.0:
        return np.zeros_like(vec)
    return vec / nrm


def _polyline_lengths(points: np.ndarray) -> np.ndarray:
    if points.shape[0] < 2:
        return np.zeros((0,), dtype=float)
    return np.linalg.norm(np.diff(points, axis=0), axis=1)


def _is_closed(points: np.ndarray, tol: float) -> bool:
    if points.shape[0] < 3:
        return False
    return float(np.linalg.norm(points[0] - points[-1])) <= tol


def _polygon_area(points: np.ndarray) -> float:
    if points.shape[0] < 3:
        return 0.0
    xy = points
    if _is_closed(xy, tol=1.0e-12):
        xy = xy[:-1]
    x = xy[:, 0]
    y = xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _split_by_mask(points: np.ndarray, mask: np.ndarray) -> list[np.ndarray]:
    segments: list[np.ndarray] = []
    if points.shape[0] == 0:
        return segments

    start = None
    for i, keep in enumerate(mask):
        if keep and start is None:
            start = i
        if (not keep or i == mask.size - 1) and start is not None:
            end = i if not keep else i + 1
            seg = points[start:end]
            if seg.shape[0] >= 2:
                segments.append(seg)
            start = None
    return segments


def _stitch_segments_rz(segments: list[np.ndarray], tol: float) -> list[np.ndarray]:
    remaining = [seg.copy() for seg in segments if seg.shape[0] >= 2]
    stitched: list[np.ndarray] = []

    def endpoints(seg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return seg[0], seg[-1]

    while remaining:
        cur = remaining.pop()
        changed = True
        while changed and remaining:
            changed = False
            cur_start, cur_end = endpoints(cur)

            best_k = None
            best_mode = None
            best_dist = None
            for k, seg in enumerate(remaining):
                s_start, s_end = endpoints(seg)
                candidates = [
                    ("append", float(np.linalg.norm(cur_end - s_start))),
                    ("append_rev", float(np.linalg.norm(cur_end - s_end))),
                    ("prepend", float(np.linalg.norm(cur_start - s_end))),
                    ("prepend_rev", float(np.linalg.norm(cur_start - s_start))),
                ]
                mode, dist = min(candidates, key=lambda x: x[1])
                if dist <= tol and (best_dist is None or dist < best_dist):
                    best_k = k
                    best_mode = mode
                    best_dist = dist

            if best_k is None or best_mode is None:
                continue

            seg = remaining.pop(best_k)
            if best_mode == "append":
                cur = np.vstack((cur, seg[1:]))
            elif best_mode == "append_rev":
                cur = np.vstack((cur, seg[::-1][1:]))
            elif best_mode == "prepend":
                cur = np.vstack((seg[:-1], cur))
            elif best_mode == "prepend_rev":
                cur = np.vstack((seg[::-1][:-1], cur))
            else:
                raise RuntimeError("unreachable")
            changed = True

        stitched.append(cur)

    return stitched


def slice_mesh_rz(
    mesh: object,
    phi: float,
    *,
    axis_xy: tuple[float, float] = (0.0, 0.0),
    stitch_tol: float = 1.0e-6,
) -> list[np.ndarray]:
    """
    Slice a triangle mesh with a constant-toroidal-angle plane.

    Returns RZ polylines in the half-plane corresponding to the given phi.

    Args:
        mesh: A trimesh.Trimesh object
        phi: Toroidal angle in radians
        axis_xy: Machine axis location in Cartesian XY
        stitch_tol: Tolerance for stitching segments

    Returns:
        List of (N, 2) arrays with R, Z coordinates
    """
    import trimesh

    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("mesh must be a trimesh.Trimesh")

    axis_origin = np.array([axis_xy[0], axis_xy[1], 0.0], dtype=float)
    e_r = np.array([np.cos(phi), np.sin(phi), 0.0], dtype=float)
    plane_normal = np.array([-np.sin(phi), np.cos(phi), 0.0], dtype=float)
    plane_origin = axis_origin

    section = mesh.section(plane_origin=plane_origin, plane_normal=plane_normal)
    if section is None:
        return []

    segments_rz: list[np.ndarray] = []
    polylines_3d: list[np.ndarray] = []
    try:
        polylines_3d = list(section.discrete)
    except (ModuleNotFoundError, Exception):
        polylines_3d = []
    if not polylines_3d:
        for ent in section.entities:
            pts = getattr(ent, "points", None)
            if pts is None:
                continue
            idx = np.asarray(pts, dtype=int)
            if idx.size < 2:
                continue
            polylines_3d.append(section.vertices[idx])

    for poly3 in polylines_3d:
        poly3_rel = poly3 - axis_origin
        r_signed = poly3_rel @ e_r
        z = poly3_rel[:, 2]
        keep = r_signed >= 0.0
        segs = _split_by_mask(np.column_stack((r_signed, z)), keep)
        segments_rz.extend(segs)

    if not segments_rz:
        return []

    return _stitch_segments_rz(segments_rz, tol=stitch_tol)


def classify_rz_contours(
    contours: list[np.ndarray], *, closed_tol: float = 1.0e-6
) -> tuple[np.ndarray, list[np.ndarray]]:
    """
    Pick an outer contour and return remaining closed contours as holes.

    If multiple closed contours exist, the outer contour is the one with the
    largest absolute signed area.
    """
    if not contours:
        raise ValueError("no contours to classify")

    closed: list[np.ndarray] = []
    open_: list[np.ndarray] = []
    for c in contours:
        if _is_closed(c, tol=closed_tol):
            closed.append(c)
        else:
            open_.append(c)

    if closed:
        areas = [abs(_polygon_area(c)) for c in closed]
        i_outer = int(np.argmax(areas))
        outer = closed[i_outer]
        holes = [c for i, c in enumerate(closed) if i != i_outer]
        return outer, holes

    outer = max(open_, key=lambda c: c.shape[0])
    return outer, []


def fill_open_contour(contour: np.ndarray, *, closed_tol: float = 1.0e-10) -> np.ndarray:
    """
    Close an open RZ polyline by adding a C1 bridge using a cubic Bezier curve.
    """
    if contour.ndim != 2 or contour.shape[1] != 2:
        raise ValueError("contour must have shape (n, 2)")
    if contour.shape[0] < 3:
        raise ValueError("contour must have at least 3 points")

    if _is_closed(contour, tol=closed_tol):
        if float(np.linalg.norm(contour[0] - contour[-1])) == 0.0:
            return contour
        return np.vstack((contour, contour[0]))

    start = contour[0]
    end = contour[-1]

    t_start = _unit(contour[1] - contour[0])
    t_end = _unit(contour[-1] - contour[-2])

    gap = float(np.linalg.norm(start - end))
    seg_lengths = _polyline_lengths(contour)
    ds = float(np.median(seg_lengths)) if seg_lengths.size else gap / 16.0
    ds = max(ds, gap / 128.0)
    n_bridge = int(max(16, np.ceil(gap / max(ds, 1.0e-12))))

    alpha = 0.35
    c1 = end + alpha * gap * t_end
    c2 = start - alpha * gap * t_start

    u = np.linspace(0.0, 1.0, n_bridge + 1)
    u = u[1:-1]
    if u.size:
        uu = u[:, None]
        bridge = (
            (1.0 - uu) ** 3 * end
            + 3.0 * (1.0 - uu) ** 2 * uu * c1
            + 3.0 * (1.0 - uu) * uu**2 * c2
            + uu**3 * start
        )
        closed = np.vstack((contour, bridge, start))
    else:
        closed = np.vstack((contour, start))
    return closed


def resample_closed_contour(contour: np.ndarray, n: int) -> np.ndarray:
    """
    Resample a closed contour to n evenly-spaced points.
    """
    if n < 4:
        raise ValueError("n must be >= 4")
    if not _is_closed(contour, tol=1.0e-8):
        raise ValueError("contour must be closed for resampling")

    pts = contour
    if float(np.linalg.norm(pts[0] - pts[-1])) == 0.0:
        pts = pts[:-1]
    ds = _polyline_lengths(np.vstack((pts, pts[0])))
    s = np.concatenate(([0.0], np.cumsum(ds)))
    total = float(s[-1])
    if total == 0.0:
        raise ValueError("degenerate contour")

    s_targets = np.linspace(0.0, total, n, endpoint=False)
    x = np.interp(s_targets, s, np.r_[pts[:, 0], pts[0, 0]])
    y = np.interp(s_targets, s, np.r_[pts[:, 1], pts[0, 1]])
    resampled = np.column_stack((x, y))
    return np.vstack((resampled, resampled[0]))


@dataclass(frozen=True)
class BoundarySlice:
    """A single toroidal slice of a wall boundary."""

    axis_xy: tuple[float, float]
    phi: float
    outer_filled: np.ndarray
    outer_original_on_filled: np.ndarray
    hole_spans: list[tuple[int, int]]
    hole_contours: list[np.ndarray]


def extract_boundary_slices(
    stl_path: Path | str,
    *,
    n_phi: int = 32,
    n_boundary_points: int = 512,
    axis_xy: tuple[float, float] | None = None,
    stitch_tol: float = 1.0e-6,
    phi_vals: np.ndarray | None = None,
) -> list[BoundarySlice]:
    """
    Extract RZ boundary slices from an STL mesh at multiple toroidal angles.

    Args:
        stl_path: Path to the STL file
        n_phi: Number of toroidal slices
        n_boundary_points: Points per boundary contour
        axis_xy: Machine axis location (auto-detected if None)
        stitch_tol: Tolerance for segment stitching
        phi_vals: Explicit toroidal angles (overrides n_phi)

    Returns:
        List of BoundarySlice objects
    """
    import trimesh

    stl_path = Path(stl_path)
    mesh = trimesh.load_mesh(stl_path, process=False)
    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("expected a single STL mesh")

    if axis_xy is None:
        axis_xy = (
            float(np.mean(mesh.vertices[:, 0])),
            float(np.mean(mesh.vertices[:, 1])),
        )

    if phi_vals is None:
        phi_vals = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    else:
        phi_vals = np.asarray(phi_vals, dtype=float)
        if phi_vals.ndim != 1 or phi_vals.size < 2:
            raise ValueError("phi_vals must be a 1D array with at least 2 values")
        n_phi = int(phi_vals.size)

    slices: list[BoundarySlice] = []
    for phi in phi_vals:
        phi_used = float(phi)
        contours = slice_mesh_rz(mesh, phi_used, axis_xy=axis_xy, stitch_tol=stitch_tol)
        if not contours:
            phi_used = float((phi_used + np.pi) % (2.0 * np.pi))
            contours = slice_mesh_rz(mesh, phi_used, axis_xy=axis_xy, stitch_tol=stitch_tol)
        if not contours:
            raise RuntimeError(f"no section contours found at phi={phi}")

        outer, holes = classify_rz_contours(contours)
        outer_filled = fill_open_contour(outer)
        outer_filled = resample_closed_contour(outer_filled, n_boundary_points)

        if _is_closed(outer, tol=1.0e-6):
            outer_on_filled = outer_filled.copy()
            hole_spans: list[tuple[int, int]] = []
        else:
            outer_on_filled, hole_spans = _nan_gaps_from_original(outer_filled, outer)

        slices.append(
            BoundarySlice(
                axis_xy=axis_xy,
                phi=phi_used,
                outer_filled=outer_filled,
                outer_original_on_filled=outer_on_filled,
                hole_spans=hole_spans,
                hole_contours=holes,
            )
        )

    return slices


def _nan_gaps_from_original(
    filled: np.ndarray, original: np.ndarray, *, dist_factor: float = 3.0
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Mark regions in filled that are far from original with NaN."""
    from scipy.spatial import cKDTree

    pts = filled[:-1]
    tree = cKDTree(original)

    seg_lengths = _polyline_lengths(original)
    ds = float(np.median(seg_lengths)) if seg_lengths.size else 0.0
    ds = max(ds, 1.0e-12)
    thresh = dist_factor * ds

    dist, _ = tree.query(pts, k=1)
    keep = dist <= thresh

    out = filled.copy()
    out[:-1][~keep] = np.nan

    holes: list[tuple[int, int]] = []
    in_gap = False
    start = 0
    for i, ok in enumerate(keep, start=1):
        if (not ok) and (not in_gap):
            in_gap = True
            start = i
        if ok and in_gap:
            holes.append((start, i - 1))
            in_gap = False
    if in_gap:
        holes.append((start, keep.size))

    return out, holes


def curve_from_contour(contour: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """
    Create a parametric curve function from a closed contour.

    Args:
        contour: Closed (N, 2) array of R, Z points

    Returns:
        Function mapping t in [0, 2*pi) to (R, Z)
    """
    pts = contour
    if float(np.linalg.norm(pts[0] - pts[-1])) == 0.0:
        pts = pts[:-1]
    ds = _polyline_lengths(np.vstack((pts, pts[0])))
    s = np.concatenate(([0.0], np.cumsum(ds)))
    total = float(s[-1])
    if total == 0.0:
        raise ValueError("degenerate contour")

    def curve(t: np.ndarray) -> np.ndarray:
        tt = np.asarray(t, dtype=float)
        tt = np.mod(tt, 2.0 * np.pi)
        s_target = tt / (2.0 * np.pi) * total
        r = np.interp(s_target, s, np.r_[pts[:, 0], pts[0, 0]])
        z = np.interp(s_target, s, np.r_[pts[:, 1], pts[0, 1]])
        return np.array([r, z])

    return curve


def export_mesh_to_stl(mesh: object, path: Path | str) -> Path:
    """Export a trimesh mesh to STL file."""
    import trimesh

    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("mesh must be a trimesh.Trimesh")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(str(path))
    return path


def convert_stl_units(
    input_path: Path | str,
    output_path: Path | str,
    *,
    from_units: str,
    to_units: str,
) -> Path:
    """
    Convert an STL file between mm and m.

    STL has no unit metadata; this applies a pure geometric scale.
    """
    import trimesh

    inp = Path(input_path)
    out = Path(output_path)

    from_u = from_units.strip().lower()
    to_u = to_units.strip().lower()

    if from_u not in {"m", "mm"}:
        raise ValueError("from_units must be m or mm")
    if to_u not in {"m", "mm"}:
        raise ValueError("to_units must be m or mm")

    if from_u == to_u:
        scale = 1.0
    elif from_u == "mm" and to_u == "m":
        scale = 1.0e-3
    elif from_u == "m" and to_u == "mm":
        scale = 1.0e3
    else:
        raise RuntimeError("unreachable units conversion")

    mesh = trimesh.load_mesh(inp, process=False)
    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError("expected a single STL mesh")

    if scale != 1.0:
        mesh = mesh.copy()
        mesh.apply_scale(scale)

    out.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(str(out))
    return out
