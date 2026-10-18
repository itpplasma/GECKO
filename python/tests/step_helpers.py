"""ALPES-style rectangular lofts and an independent finite-segment field oracle.

Generator adapted from ALPES CadQueryMultiCoilGenerator.py and
windingCoordinateGenerator.py (MIT), revision
 eb31dc8f4725a5bc355a558b4357f9234c312264. See data/README.md.
"""

import numpy as np


def loft_coil(points, width=40.0, height=67.0, *, fuse=True):
    """Reproduce ALPES customOrientePlanes and sequential loft/union logic."""
    import cadquery as cq

    tangent = np.roll(points, -1, axis=0) - np.roll(points, 1, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    # ALPES uses Y as the steering direction for coil 0.
    delta = points.mean(axis=0) - points
    radial = delta - delta[:, 1, None] * np.array([0.0, 1.0, 0.0])
    radial /= np.linalg.norm(radial, axis=1)[:, None]
    xdir = np.cross(tangent, radial)
    xdir /= np.linalg.norm(xdir, axis=1)[:, None]
    profiles = []
    for point, x, normal in zip(points, xdir, tangent):
        plane = cq.Plane(origin=tuple(point), xDir=tuple(x), normal=tuple(normal))
        profiles.append(cq.Workplane(plane).rect(width, height))
    # CadQuery booleans default to parallel execution. Bound the OCCT pool.
    from OCP.OSD import OSD_Parallel, OSD_ThreadPool

    OSD_Parallel.SetUseOcctThreads_s(True)
    pool = OSD_ThreadPool.DefaultPool_s()
    pool.Init(8)
    pool.SetNbDefaultThreadsToLaunch(8)
    result = None
    sides = []
    for i in range(len(points)):
        wp = cq.Workplane()
        wp.ctx.pendingWires.extend([profiles[i - 1].val(), profiles[i].val()])
        segment = wp.loft(combine=True)
        if fuse:
            result = segment if result is None else result.union(segment)
        else:
            for face in segment.val().Faces():
                center = np.array(face.Center().toTuple())
                if (
                    min(
                        np.linalg.norm(center - points[i - 1]),
                        np.linalg.norm(center - points[i]),
                    )
                    > 1e-6
                ):
                    sides.append(face)
    return result.val() if fuse else cq.Solid.makeSolid(cq.Shell.makeShell(sides))


def field(points, line, current):
    """Exact Biot-Savart integral over straight segments, in SI units."""
    a = line - points[:, None, :]
    b = np.roll(line, -1, axis=0) - points[:, None, :]
    ra = np.linalg.norm(a, axis=2)
    rb = np.linalg.norm(b, axis=2)
    product = ra * rb
    factor = (ra + rb) / (product * (product + np.sum(a * b, axis=2)))
    return 1e-7 * current * np.sum(np.cross(a, b) * factor[:, :, None], axis=1)


def nearest_on_polyline(points, line):
    """Return closest projections on all segments (including closure)."""
    segments = np.roll(line, -1, axis=0) - line
    delta = points[:, None, :] - line
    t = np.sum(delta * segments, axis=2) / np.sum(segments**2, axis=1)
    projected = line + np.clip(t, 0, 1)[:, :, None] * segments
    distance = np.linalg.norm(points[:, None, :] - projected, axis=2)
    index = np.argmin(distance, axis=1)
    return projected[np.arange(len(points)), index], distance.min(axis=1)
