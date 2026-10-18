"""STEP winding packs to filaments, with CadQuery loaded only on demand.

Supports closed rectangular section lofts with a conforming quadrilateral
boundary grid (including the ALPES generator). No CAD current direction is
assumed: positive current follows a canonical geometric orientation.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from .coils import CoilsFile, Filament


class StepCoilError(ValueError):
    """Invalid, open, self-intersecting or unsupported winding-pack geometry."""


def _section_centers(solid) -> np.ndarray:
    """Separate the two edge families of a closed quadrilateral tube grid."""
    import cadquery as cq

    edges = solid.Edges()
    vertices = solid.Vertices()
    edge_id = {edge: i for i, edge in enumerate(edges)}
    vertex_id = {vertex: i for i, vertex in enumerate(vertices)}
    ends = [{vertex_id[v] for v in edge.Vertices()} for edge in edges]
    if any(len(pair) != 2 for pair in ends):
        raise StepCoilError("Unsupported topology: expected segmented loft edges")
    incident = [set() for _ in vertices]
    adjacent = [set() for _ in edges]
    faces = []
    parent = list(range(len(edges)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def join(i, j):
        parent[root(i)] = root(j)

    for i, pair in enumerate(ends):
        for v in pair:
            incident[v].add(i)
    for face in solid.Faces():
        ids = [edge_id[edge] for edge in face.Edges()]
        if len(ids) != 4 or len(face.Wires()) != 1:
            raise StepCoilError(
                "Unsupported topology: expected quadrilateral loft faces"
            )
        faces.append(ids)
        for i in ids:
            for j in ids:
                if i != j:
                    if ends[i] & ends[j]:
                        adjacent[i].add(j)
                    else:
                        join(i, j)
    for ids in incident:
        if len(ids) != 4:
            raise StepCoilError(
                "Winding pack is not a closed loft grid (open or capped)"
            )
        for i in ids:
            for j in ids - adjacent[i] - {i}:
                join(i, j)
    families = {}
    for i in range(len(edges)):
        families.setdefault(root(i), set()).add(i)
    candidates = []
    for family in families.values():
        remaining = set(family)
        loops = []
        membership = {}
        while remaining:
            pending = [min(remaining)]
            component = set()
            while pending:
                i = pending.pop()
                if i in component:
                    continue
                component.add(i)
                for v in ends[i]:
                    pending.extend((incident[v] & family) - component)
            remaining -= component
            if any(len(incident[v] & family) != 2 for i in component for v in ends[i]):
                break
            for i in component:
                membership[i] = len(loops)
            loops.append(component)
        else:
            if len(loops) < 3 or len(membership) != len(family):
                continue
            # Sections must be planar rectangles, not longitudinal rails.
            centers = []
            try:
                for loop in loops:
                    if len(loop) != 4:
                        raise ValueError("Not a rectangular section")
                    wire = cq.Wire.assembleEdges([edges[i] for i in sorted(loop)])
                    section = cq.Face.makeFromWires(wire)
                    if (
                        not wire.IsClosed()
                        or not section.isValid()
                        or section.geomType() != "PLANE"
                    ):
                        raise ValueError("Invalid section")
                    corners = np.array([v.toTuple() for v in wire.Vertices()])
                    center = corners.mean(axis=0)
                    # A planar parallelogram with orthogonal adjacent sides.
                    radii = np.linalg.norm(corners - center, axis=1)
                    if np.ptp(radii) > 1e-6 * radii.max():
                        raise ValueError("Section is not rectangular")
                    centers.append(section.Center().toTuple())
            except (ValueError, RuntimeError):
                continue
            neighbors = [set() for _ in loops]
            for ids in faces:
                pair = {membership[i] for i in ids if i in family}
                if len(pair) != 2:
                    break
                a, b = pair
                neighbors[a].add(b)
                neighbors[b].add(a)
            else:
                if any(len(row) != 2 for row in neighbors):
                    continue
                order = [0]
                previous, current = -1, 0
                while True:
                    following = min(neighbors[current] - {previous})
                    if following == 0:
                        break
                    if following in order:
                        break
                    order.append(following)
                    previous, current = current, following
                if len(order) == len(loops):
                    perimeter = sum(edges[i].Length() for i in family) / len(loops)
                    candidates.append((perimeter, np.array(centers)[order]))
    if not candidates:
        raise StepCoilError(
            "Unsupported or non-closed winding pack: need a closed chain of "
            "planar rectangular sections joined by quadrilateral loft faces"
        )
    return min(candidates, key=lambda item: item[0])[1]


def _resample(points: np.ndarray, n_points: int) -> np.ndarray:
    # Fix the sign of the largest projected oriented area, then the start point.
    area = np.cross(points, np.roll(points, -1, axis=0)).sum(axis=0)
    if np.linalg.norm(area) < 1e-12:
        raise StepCoilError("Degenerate centerline: orientation is undefined")
    if area[np.argmax(np.abs(area))] < 0:
        points = points[::-1]
    start = np.lexsort(points.T[::-1])[0]
    points = np.roll(points, -start, axis=0)
    closed = np.vstack([points, points[0]])
    lengths = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    if np.any(lengths <= 1e-12):
        raise StepCoilError("Degenerate centerline: coincident section centers")
    arc = np.r_[0.0, np.cumsum(lengths)]
    sample = np.linspace(0.0, arc[-1], n_points + 1)
    result = np.column_stack([np.interp(sample, arc, col) for col in closed.T])
    result[-1] = result[0]
    return result


def extract_step_centerlines(
    filename: str | Path, *, n_points: int = 256, units: str = "mm"
) -> list[np.ndarray]:
    """Extract one closed centerline per STEP solid, returned in meters.

    ``n_points`` is the number of unique, uniform arc-length samples; the first
    point is repeated at the end. ``units`` selects the OpenCASCADE working unit
    (mm or m), respecting the STEP header, rather than overriding its units.
    Solids retain STEP traversal order. Unsupported geometry is rejected rather
    than approximated silently. Positive orientation has the largest component
    of the polygon's oriented area vector positive; the start is lexicographic.
    """
    if units not in ("mm", "m"):
        raise ValueError("units must be 'mm' or 'm'")
    if isinstance(n_points, bool) or not isinstance(n_points, int) or n_points < 8:
        raise ValueError("n_points must be an integer >= 8")
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        import cadquery as cq
    except ImportError as exc:
        raise ImportError("STEP import requires gecko[step] (CadQuery)") from exc
    try:
        imported = cq.importers.importStep(str(path), unit=units.upper())
        shape = cq.Compound.makeCompound(imported.vals())
    except (ValueError, RuntimeError) as exc:
        raise StepCoilError(f"Cannot read STEP file: {path.name}") from exc
    solids = shape.Solids()
    if not solids or sum(len(s.Faces()) for s in solids) != len(shape.Faces()):
        raise StepCoilError(
            "STEP must contain closed solids, with no loose faces/shells"
        )
    from OCP.BOPAlgo import BOPAlgo_ArgumentAnalyzer

    lines = []
    for index, solid in enumerate(solids, 1):
        try:
            if not solid.isValid() or any(not s.Closed() for s in solid.Shells()):
                raise StepCoilError("Invalid, non-closed or self-intersecting solid")
            analyzer = BOPAlgo_ArgumentAnalyzer()
            analyzer.SetShape1(solid.wrapped)
            analyzer.SelfInterMode = True
            analyzer.SetRunParallel(False)
            analyzer.SetUseOBB(True)
            analyzer.Perform()
            if analyzer.HasFaulty():
                raise StepCoilError("Self-intersecting winding-pack solid")
            points = _section_centers(solid) * (1e-3 if units == "mm" else 1.0)
            lines.append(_resample(points, n_points))
        except StepCoilError as exc:
            raise StepCoilError(f"Coil {index}: {exc}") from exc
    return lines


def coils_from_step(
    filename: str | Path,
    *,
    currents: float | Sequence[float],
    n_points: int = 256,
    units: str = "mm",
    reverse: Sequence[int] = (),
) -> CoilsFile:
    """Build filaments using signed currents in amperes (scalar or per solid).

    ``reverse`` contains zero-based solid indices whose canonical orientation
    should be reversed. Geometry alone cannot determine physical current sense.
    No symmetry copies are generated; periods=1 and mirror=NIL.
    """
    lines = extract_step_centerlines(filename, n_points=n_points, units=units)
    values = np.asarray(currents, dtype=float)
    if values.ndim == 0:
        values = np.full(len(lines), float(values))
    if values.shape != (len(lines),) or not np.isfinite(values).all():
        raise ValueError("Supply one finite current per solid, or a scalar current")
    if any(
        isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(lines)
        for i in reverse
    ):
        raise ValueError("reverse must contain zero-based solid indices")
    filaments = []
    for i, (line, current) in enumerate(zip(lines, values)):
        if i in reverse:
            line = line[::-1].copy()
        filaments.append(Filament(line, np.full(len(line), current)))
    return CoilsFile(periods=1, mirror="NIL", filaments=filaments)


def write_step_coils(
    filename: str | Path,
    output: str | Path,
    *,
    currents: float | Sequence[float],
    n_points: int = 256,
    units: str = "mm",
    reverse: Sequence[int] = (),
) -> CoilsFile:
    """Extract and write MAKEGRID coils.* in meters and amperes.

    A repeated closing coordinate carries zero current and the coil group tag,
    as required by MAKEGRID. Zero-current coils are omitted from the file.
    The returned CoilsFile retains all solids, including zero-current ones.
    """
    coils = coils_from_step(
        filename, currents=currents, n_points=n_points, units=units, reverse=reverse
    )
    if not any(fil.current[0] != 0 for fil in coils.filaments):
        raise ValueError("At least one nonzero current is required for coils.* output")
    with Path(output).open("w", encoding="utf-8") as stream:
        stream.write("periods 1\nbegin filament\nmirror NIL\n")
        for group, fil in enumerate(coils.filaments, 1):
            current = float(fil.current[0])
            if current == 0:
                continue
            stream.writelines(
                f"{x:.14E} {y:.14E} {z:.14E} {current:.14E}\n"
                for x, y, z in fil.coords[:-1]
            )
            x, y, z = fil.coords[-1]
            stream.write(f"{x:.14E} {y:.14E} {z:.14E} 0 {group} coil_{group}\n")
        stream.write("end\n")
    return coils
