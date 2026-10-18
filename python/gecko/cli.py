"""
Command-line interface for GECKO.

Usage:
    gecko view wall.stl coils.scaled
    gecko wall-from-vmec wout.nc --offset 0.25 -o wall.stl
    gecko coils-info coils.scaled
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def cmd_view(args: argparse.Namespace) -> int:
    """View geometry files interactively."""
    from .viewer import quick_view

    cross_section = None
    if args.coil_width and args.coil_height:
        cross_section = (args.coil_width, args.coil_height)

    if args.stl_scale is not None:
        stl_scale = float(args.stl_scale)
    else:
        stl_units = str(args.stl_units).strip().lower()
        if stl_units == "m":
            stl_scale = 1.0
        elif stl_units == "mm":
            stl_scale = 1.0e-3
        else:
            raise ValueError("stl_units must be m or mm")

    quick_view(
        *args.files,
        coil_cross_section=cross_section,
        show=True,
        stl_scale=stl_scale,
    )
    return 0


def cmd_wall_from_vmec(args: argparse.Namespace) -> int:
    """Generate wall STL from VMEC wout file."""
    from .stl import export_mesh_to_stl
    from .wall import offset_wall_mesh, wall_mesh_from_wout

    if args.offset > 0:
        mesh = offset_wall_mesh(
            args.wout,
            args.offset,
            n_theta=args.n_theta,
            n_zeta=args.n_zeta,
        )
    else:
        mesh = wall_mesh_from_wout(
            args.wout,
            n_theta=args.n_theta,
            n_zeta=args.n_zeta,
        )

    output = args.output or Path(args.wout).with_suffix(".stl")
    export_mesh_to_stl(mesh, output)
    print(f"Wrote {output}")
    return 0


def cmd_coils_info(args: argparse.Namespace) -> int:
    """Show information about a coils file."""
    from libneo.coils import CoilsFile

    coils = CoilsFile.from_file(str(args.coils))
    print(f"File: {args.coils}")
    print(f"Periods: {coils.periods}")
    print(f"Mirror: {coils.mirror}")
    print(f"Number of coils: {len(coils.filaments)}")
    print(f"Total points: {sum(len(f.coords) for f in coils.filaments)}")

    all_pts = []
    for fil in coils.filaments:
        all_pts.append(fil.coords)
    if all_pts:
        import numpy as np
        stacked = np.vstack(all_pts)
        bbox_min = stacked.min(axis=0)
        bbox_max = stacked.max(axis=0)
        print(f"Bounding box:")
        print(f"  Min: ({bbox_min[0]:.3f}, {bbox_min[1]:.3f}, {bbox_min[2]:.3f})")
        print(f"  Max: ({bbox_max[0]:.3f}, {bbox_max[1]:.3f}, {bbox_max[2]:.3f})")

    for i, fil in enumerate(coils.filaments):
        cur = float(fil.current.mean())
        print(f"  Coil {i+1}: {len(fil.coords)} points, I = {cur:.3e} A")

    return 0


def cmd_coils_plot(args: argparse.Namespace) -> int:
    """Plot coils in 3D."""
    from libneo.coils import CoilsFile

    coils = CoilsFile.from_file(str(args.coils))
    fig, _ = coils.plot(show=not args.output, show_legend=True)

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"Wrote {args.output}")

    return 0


def cmd_stl_convert(args: argparse.Namespace) -> int:
    """Convert STL units between mm and m."""
    from .stl import convert_stl_units

    out = convert_stl_units(
        args.input,
        args.output,
        from_units=args.from_units,
        to_units=args.to_units,
    )
    print(f"Wrote {out}")
    return 0


def cmd_view_results(args: argparse.Namespace) -> int:
    """View SIMPLE results.nc wall-hit points overlaid on an STL wall."""
    from .simple_results import load_simple_wall_hits
    from .viewer import Viewer

    if args.stl_scale is not None:
        stl_scale = float(args.stl_scale)
    else:
        stl_units = str(args.stl_units).strip().lower()
        if stl_units == "m":
            stl_scale = 1.0
        elif stl_units == "mm":
            stl_scale = 1.0e-3
        else:
            raise ValueError("stl_units must be m or mm")

    hits = load_simple_wall_hits(args.results)

    viewer = Viewer(off_screen=bool(args.off_screen))
    viewer.add_axes()

    if args.wall is not None:
        viewer.add_wall(args.wall, opacity=float(args.wall_opacity), scale=stl_scale)

    if args.coils is not None:
        viewer.add_coils(
            args.coils,
            cross_section=(args.coil_width, args.coil_height)
            if args.show_coil_mesh
            else None,
            show_centerline=bool(args.show_coil_centerline),
        )

    if args.color_by == "cos_incidence" and hits.cos_incidence is not None:
        viewer.add_points(
            hits.points_m,
            scalars=hits.cos_incidence,
            cmap=str(args.cmap),
            point_size=float(args.point_size),
        )
        viewer.add_scalar_bar(title="cos_incidence")
    elif args.color_by == "angle_rad" and hits.angle_rad is not None:
        viewer.add_points(
            hits.points_m,
            scalars=hits.angle_rad,
            cmap=str(args.cmap),
            point_size=float(args.point_size),
        )
        viewer.add_scalar_bar(title="angle_rad")
    else:
        viewer.add_points(
            hits.points_m,
            color=str(args.point_color),
            point_size=float(args.point_size),
        )

    if args.screenshot is not None:
        viewer.screenshot(args.screenshot)
    if not args.off_screen:
        viewer.show()
    return 0


def cmd_coils_from_step(args: argparse.Namespace) -> int:
    """Extract STEP winding packs and write filaments in SI units."""
    from .step import StepCoilError, write_step_coils

    currents = args.currents[0] if len(args.currents) == 1 else args.currents
    try:
        coils = write_step_coils(
            args.input, args.output, currents=currents, n_points=args.n_points,
            units=args.units, reverse=args.reverse,
        )
    except (StepCoilError, ValueError, OSError, ImportError) as exc:
        print(f"STEP import failed: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {args.output}: {len(coils.filaments)} coils (m, A)")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="gecko",
        description="GECKO - Geometry, Engineering, and CAD for Kinetic modeling",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # view command
    p_view = subparsers.add_parser("view", help="View geometry files interactively")
    p_view.add_argument("files", nargs="+", type=Path, help="Files to view")
    p_view.add_argument("--coil-width", type=float, default=0.4)
    p_view.add_argument("--coil-height", type=float, default=0.8)
    p_view.add_argument(
        "--stl-units",
        choices=["m", "mm"],
        default="m",
        help="Interpret STL coordinates in these units (no auto-detect).",
    )
    p_view.add_argument(
        "--stl-scale",
        type=float,
        default=None,
        help="Additional STL scale factor (overrides --stl-units).",
    )
    p_view.set_defaults(func=cmd_view)

    # wall-from-vmec command
    p_wall = subparsers.add_parser("wall-from-vmec", help="Generate wall STL from VMEC")
    p_wall.add_argument("wout", type=Path, help="VMEC wout file")
    p_wall.add_argument("-o", "--output", type=Path, help="Output STL file")
    p_wall.add_argument("--offset", type=float, default=0.0, help="Offset from LCFS")
    p_wall.add_argument("--n-theta", type=int, default=128)
    p_wall.add_argument("--n-zeta", type=int, default=96)
    p_wall.set_defaults(func=cmd_wall_from_vmec)

    p_step = subparsers.add_parser(
        "coils-from-step", help="Extract filament centerlines from STEP winding packs"
    )
    p_step.add_argument("input", type=Path)
    p_step.add_argument("-o", "--output", type=Path, required=True)
    p_step.add_argument(
        "--currents", type=float, nargs="+", required=True,
        help="Signed amperes: one value for all solids, or one per STEP solid",
    )
    p_step.add_argument("--n-points", type=int, default=256)
    p_step.add_argument(
        "--units", choices=["mm", "m"], default="mm",
        help="CAD working unit; STEP header units are respected; output is meters",
    )
    p_step.add_argument(
        "--reverse", type=int, nargs="*", default=[],
        help="Zero-based solid indices to reverse from canonical orientation",
    )
    p_step.set_defaults(func=cmd_coils_from_step)

    # coils-info command
    p_info = subparsers.add_parser("coils-info", help="Show coils file information")
    p_info.add_argument("coils", type=Path, help="Coils file")
    p_info.set_defaults(func=cmd_coils_info)

    # coils-plot command
    p_plot = subparsers.add_parser("coils-plot", help="Plot coils in 3D")
    p_plot.add_argument("coils", type=Path, help="Coils file")
    p_plot.add_argument("-o", "--output", type=Path, help="Output image file")
    p_plot.set_defaults(func=cmd_coils_plot)

    p_stl = subparsers.add_parser("stl-convert", help="Convert STL units (mm <-> m)")
    p_stl.add_argument("input", type=Path, help="Input STL file")
    p_stl.add_argument("output", type=Path, help="Output STL file")
    p_stl.add_argument("--from-units", choices=["m", "mm"], required=True)
    p_stl.add_argument("--to-units", choices=["m", "mm"], required=True)
    p_stl.set_defaults(func=cmd_stl_convert)

    p_results = subparsers.add_parser(
        "view-results",
        help="View SIMPLE results.nc (wall_hit points) overlaid on an STL wall",
    )
    p_results.add_argument("results", type=Path, help="SIMPLE results.nc file")
    p_results.add_argument("--wall", type=Path, default=None, help="Wall STL file")
    p_results.add_argument("--coils", type=Path, default=None, help="Optional coils file")
    p_results.add_argument("--wall-opacity", type=float, default=0.6)
    p_results.add_argument("--point-size", type=float, default=6.0)
    p_results.add_argument("--point-color", type=str, default="red")
    p_results.add_argument("--cmap", type=str, default="viridis")
    p_results.add_argument(
        "--color-by",
        choices=["none", "cos_incidence", "angle_rad"],
        default="none",
        help="Color hit points by a scalar from results.nc (if present).",
    )
    p_results.add_argument(
        "--stl-units",
        choices=["m", "mm"],
        default="m",
        help="Interpret STL coordinates in these units (no auto-detect).",
    )
    p_results.add_argument(
        "--stl-scale",
        type=float,
        default=None,
        help="Additional STL scale factor (overrides --stl-units).",
    )
    p_results.add_argument("--off-screen", action="store_true")
    p_results.add_argument("--screenshot", type=Path, default=None)
    p_results.add_argument("--coil-width", type=float, default=0.4)
    p_results.add_argument("--coil-height", type=float, default=0.8)
    p_results.add_argument("--show-coil-mesh", action="store_true")
    p_results.add_argument("--show-coil-centerline", action="store_true")
    p_results.set_defaults(func=cmd_view_results)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
