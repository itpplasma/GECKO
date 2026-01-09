"""
Coil file handling for fusion device geometry.

Supports STELLOPT/MAKEGRID coil file format and provides visualization.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Filament:
    """A single coil filament with coordinates and current."""

    coords: np.ndarray  # shape (N, 3) -> columns x, y, z
    current: np.ndarray  # shape (N,)


@dataclass(frozen=True)
class CoilSet:
    """Collection of coil filaments with optional metadata."""

    periods: int | None
    mirror: str | None
    filaments: list[Filament]


def parse_coils_file(path: Path | str) -> CoilSet:
    """
    Parse a STELLOPT/MAKEGRID coils file.

    Format:
        - 'periods N' line specifies field period symmetry
        - 'mirror ON/OFF/NIL' specifies up-down symmetry
        - 'begin filament' starts a coil block
        - 'end' or 'end filament' ends a block
        - Within a block: x y z current

    Args:
        path: Path to coils file

    Returns:
        CoilSet with parsed data
    """
    path = Path(path)
    periods: int | None = None
    mirror: str | None = None
    filaments: list[Filament] = []
    in_filament = False
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    Is: list[float] = []

    with path.open("r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("!"):
                continue
            for cmt in ("#", "!"):
                cpos = line.find(cmt)
                if cpos >= 0:
                    line = line[:cpos].strip()
            if not line:
                continue

            low = line.lower()
            if low.startswith("periods"):
                toks = line.split()
                if len(toks) > 1:
                    try:
                        periods = int(toks[1])
                    except ValueError:
                        pass
                continue
            if low.startswith("mirror"):
                toks = line.split()
                mirror = toks[1] if len(toks) > 1 else None
                continue
            if low.startswith("begin filament"):
                if in_filament and len(xs) > 0:
                    coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
                    cur = np.array(Is)
                    filaments.append(Filament(coords=coords, current=cur))
                    xs.clear()
                    ys.clear()
                    zs.clear()
                    Is.clear()
                elif in_filament and len(xs) == 0:
                    raise ValueError(f"Empty filament at {path}:{lineno}")
                in_filament = True
                continue
            if low.startswith("end filament") or low == "end":
                if not in_filament:
                    continue
                if len(xs) == 0:
                    raise ValueError(f"Empty filament at {path}:{lineno}")
                coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
                cur = np.array(Is)
                filaments.append(Filament(coords=coords, current=cur))
                xs.clear()
                ys.clear()
                zs.clear()
                Is.clear()
                in_filament = False
                continue

            if in_filament:
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(
                        f"Expected 4 numeric fields (x y z I) at {path}:{lineno}"
                    )
                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
                cur = float(parts[3])
                xs.append(x)
                ys.append(y)
                zs.append(z)
                Is.append(cur)
                continue

    if in_filament:
        if len(xs) == 0:
            raise ValueError(f"Empty filament at end of file: {path}")
        coords = np.column_stack([np.array(xs), np.array(ys), np.array(zs)])
        cur = np.array(Is)
        filaments.append(Filament(coords=coords, current=cur))

    return CoilSet(periods=periods, mirror=mirror, filaments=filaments)


class CoilsFile:
    """
    High-level interface for coil file handling.

    Supports reading, writing, plotting, and format conversion.
    """

    def __init__(
        self,
        *,
        periods: int | None,
        mirror: str | None,
        filaments: list[Filament],
    ):
        self.periods = periods
        self.mirror = mirror
        self.filaments = filaments

    @classmethod
    def from_file(cls, filename: Path | str) -> "CoilsFile":
        """Load coils from a STELLOPT/MAKEGRID format file."""
        coils = parse_coils_file(filename)
        return cls(
            periods=coils.periods, mirror=coils.mirror, filaments=coils.filaments
        )

    def write(self, filename: Path | str) -> None:
        """Write coils in STELLOPT/MAKEGRID format."""
        filename = Path(filename)
        with filename.open("w", encoding="utf-8") as f:
            if self.periods is not None:
                f.write(f"periods {self.periods}\n")
            if self.mirror is not None:
                f.write(f"mirror {self.mirror}\n")
            for fil in self.filaments:
                f.write("begin filament\n")
                for (x, y, z), I in zip(fil.coords, fil.current):
                    f.write(f"{x:.14E}   {y:.14E}   {z:.14E}   {I:.14E}\n")

    def write_simple(self, filename: Path | str) -> None:
        """
        Write coils in simple biotsavart format.

        Format:
            n_points
            x y z current
            ...
        """
        filename = Path(filename)
        with filename.open("w", encoding="utf-8") as f:
            n_points = sum(len(fil.coords) for fil in self.filaments)
            f.write(f"{n_points}\n")
            for fil in self.filaments:
                for (x, y, z), I in zip(fil.coords, fil.current):
                    f.write(f"{x:.14E}   {y:.14E}   {z:.14E}   {I:.14E}\n")

    @property
    def n_coils(self) -> int:
        """Number of coil filaments."""
        return len(self.filaments)

    @property
    def total_points(self) -> int:
        """Total number of points across all filaments."""
        return sum(len(fil.coords) for fil in self.filaments)

    def get_all_points(self) -> np.ndarray:
        """Get all coil points as a single (N, 3) array."""
        if not self.filaments:
            return np.zeros((0, 3), dtype=float)
        return np.vstack([fil.coords for fil in self.filaments])

    def bounding_box(self) -> tuple[np.ndarray, np.ndarray]:
        """Get bounding box (min, max) of all coil points."""
        pts = self.get_all_points()
        if pts.shape[0] == 0:
            return np.zeros(3), np.zeros(3)
        return pts.min(axis=0), pts.max(axis=0)

    def plot(self, *, show_legend: bool = True, show: bool = False):
        """
        Plot coils in 3D using matplotlib.

        Args:
            show_legend: Whether to show legend
            show: Whether to call plt.show()

        Returns:
            fig, ax tuple
        """
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        all_points = []
        for k, fil in enumerate(self.filaments, start=1):
            pts = fil.coords
            cur = float(fil.current.mean())
            label = f"filament_{k:02d} (I={cur:.3e} A)"
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], label=label)
            all_points.append(pts)
        if all_points:
            stacked = np.vstack(all_points)
            span = stacked.max(axis=0) - stacked.min(axis=0)
            span[span == 0.0] = 1.0
            ax.set_box_aspect(span)
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.set_zlabel("Z [m]")
        if show_legend:
            ax.legend(fontsize="small")
        plt.tight_layout()
        if show:
            plt.show()
        return fig, ax


def coil_to_tube_mesh(
    coords: np.ndarray,
    *,
    width: float,
    height: float,
    n_sides: int = 4,
) -> "trimesh.Trimesh":
    """
    Create a 3D tube mesh from coil centerline with rectangular cross-section.

    Args:
        coords: (N, 3) coil centerline points
        width: Cross-section width (horizontal)
        height: Cross-section height (vertical)
        n_sides: Number of sides (4 for rectangular)

    Returns:
        trimesh.Trimesh object
    """
    import trimesh

    n_pts = coords.shape[0]
    if n_pts < 2:
        raise ValueError("Need at least 2 points for tube mesh")

    # Compute tangent vectors
    tangent = np.zeros_like(coords)
    tangent[1:-1] = coords[2:] - coords[:-2]
    tangent[0] = coords[1] - coords[0]
    tangent[-1] = coords[-1] - coords[-2]
    tangent = tangent / (np.linalg.norm(tangent, axis=1, keepdims=True) + 1e-12)

    # Compute local coordinate frames using parallel transport
    up = np.array([0.0, 0.0, 1.0])
    vertices = []
    faces = []

    # Create cross-section profile (rectangular)
    half_w = width / 2
    half_h = height / 2
    profile = np.array([
        [-half_w, -half_h],
        [half_w, -half_h],
        [half_w, half_h],
        [-half_w, half_h],
    ])

    for i, (center, t) in enumerate(zip(coords, tangent)):
        # Compute local frame
        if abs(np.dot(t, up)) > 0.99:
            side = np.cross(t, np.array([1.0, 0.0, 0.0]))
        else:
            side = np.cross(t, up)
        side = side / (np.linalg.norm(side) + 1e-12)
        up_local = np.cross(side, t)
        up_local = up_local / (np.linalg.norm(up_local) + 1e-12)

        # Transform profile to world coordinates
        for px, py in profile:
            pt = center + px * side + py * up_local
            vertices.append(pt)

    vertices = np.array(vertices)

    # Create faces connecting adjacent profiles
    n_profile = len(profile)
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

    # Cap the ends
    def add_cap(base_idx: int, reverse: bool = False):
        center_idx = len(vertices)
        cap_verts = vertices[base_idx : base_idx + n_profile]
        center = cap_verts.mean(axis=0)
        nonlocal vertices
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
        coils: CoilsFile object
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
