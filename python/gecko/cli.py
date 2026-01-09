"""
Command-line interface for GECKO.

Usage:
    gecko view wall.stl coils.scaled
    gecko wall-from-vmec wout.nc --offset 0.25 -o wall.stl
    gecko coils-info coils.scaled
"""

from __future__ import annotations

import argparse
from pathlib import Path


def cmd_view(args: argparse.Namespace) -> int:
    """View geometry files interactively."""
    from .viewer import quick_view

    cross_section = None
    if args.coil_width and args.coil_height:
        cross_section = (args.coil_width, args.coil_height)

    quick_view(*args.files, coil_cross_section=cross_section, show=True)
    return 0


def cmd_wall_from_vmec(args: argparse.Namespace) -> int:
    """Generate wall STL from VMEC wout file."""
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
    mesh.export(str(output))
    print(f"Wrote {output}")
    return 0


def cmd_coils_info(args: argparse.Namespace) -> int:
    """Show information about a coils file."""
    from .coils import CoilsFile

    coils = CoilsFile.from_file(args.coils)
    print(f"File: {args.coils}")
    print(f"Periods: {coils.periods}")
    print(f"Mirror: {coils.mirror}")
    print(f"Number of coils: {coils.n_coils}")
    print(f"Total points: {coils.total_points}")

    bbox_min, bbox_max = coils.bounding_box()
    print(f"Bounding box:")
    print(f"  Min: ({bbox_min[0]:.3f}, {bbox_min[1]:.3f}, {bbox_min[2]:.3f})")
    print(f"  Max: ({bbox_max[0]:.3f}, {bbox_max[1]:.3f}, {bbox_max[2]:.3f})")

    for i, fil in enumerate(coils.filaments):
        cur = float(fil.current.mean())
        print(f"  Coil {i+1}: {len(fil.coords)} points, I = {cur:.3e} A")

    return 0


def cmd_coils_plot(args: argparse.Namespace) -> int:
    """Plot coils in 3D."""
    from .coils import CoilsFile

    coils = CoilsFile.from_file(args.coils)
    fig, ax = coils.plot(show=not args.output, show_legend=True)

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"Wrote {args.output}")

    return 0


def main(argv: list[str] | None = None) -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="gecko",
        description="GECKO - Geometry, Engineering, and CAD for Kinetic modeling and Optimization",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # view command
    p_view = subparsers.add_parser("view", help="View geometry files interactively")
    p_view.add_argument("files", nargs="+", type=Path, help="Files to view")
    p_view.add_argument("--coil-width", type=float, default=0.4)
    p_view.add_argument("--coil-height", type=float, default=0.8)
    p_view.set_defaults(func=cmd_view)

    # wall-from-vmec command
    p_wall = subparsers.add_parser("wall-from-vmec", help="Generate wall STL from VMEC")
    p_wall.add_argument("wout", type=Path, help="VMEC wout file")
    p_wall.add_argument("-o", "--output", type=Path, help="Output STL file")
    p_wall.add_argument("--offset", type=float, default=0.0, help="Offset from LCFS")
    p_wall.add_argument("--n-theta", type=int, default=128)
    p_wall.add_argument("--n-zeta", type=int, default=96)
    p_wall.set_defaults(func=cmd_wall_from_vmec)

    # coils-info command
    p_info = subparsers.add_parser("coils-info", help="Show coils file information")
    p_info.add_argument("coils", type=Path, help="Coils file")
    p_info.set_defaults(func=cmd_coils_info)

    # coils-plot command
    p_plot = subparsers.add_parser("coils-plot", help="Plot coils in 3D")
    p_plot.add_argument("coils", type=Path, help="Coils file")
    p_plot.add_argument("-o", "--output", type=Path, help="Output image file")
    p_plot.set_defaults(func=cmd_coils_plot)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
