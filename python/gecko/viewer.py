"""
Interactive 3D visualization for fusion geometry using PyVista.

Provides a unified viewer for walls, coils, and other geometry.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


class Viewer:
    """
    Interactive 3D viewer for fusion device geometry.

    Example:
        viewer = Viewer()
        viewer.add_wall("wall.stl")
        viewer.add_coils("coils.scaled", cross_section=(0.4, 0.8))
        viewer.show()
    """

    def __init__(self, *, background: str = "white"):
        """
        Initialize viewer.

        Args:
            background: Background color
        """
        import pyvista as pv

        self._plotter = pv.Plotter()
        self._plotter.set_background(background)
        self._meshes: list[Any] = []

    def add_wall(
        self,
        source: Path | str | Any,
        *,
        color: str = "lightgray",
        opacity: float = 0.7,
        show_edges: bool = False,
    ) -> "Viewer":
        """
        Add a wall mesh to the viewer.

        Args:
            source: STL file path or trimesh/pyvista mesh object
            color: Mesh color
            opacity: Mesh opacity (0-1)
            show_edges: Whether to show mesh edges

        Returns:
            self for chaining
        """
        import pyvista as pv

        if isinstance(source, (str, Path)):
            mesh = pv.read(str(source))
        else:
            import trimesh

            if isinstance(source, trimesh.Trimesh):
                mesh = pv.wrap(source)
            else:
                mesh = source

        self._plotter.add_mesh(
            mesh,
            color=color,
            opacity=opacity,
            show_edges=show_edges,
        )
        self._meshes.append(mesh)
        return self

    def add_coils(
        self,
        source: Path | str | Any,
        *,
        cross_section: tuple[float, float] | None = None,
        color: str = "copper",
        show_centerline: bool = False,
        centerline_color: str = "blue",
        centerline_width: float = 2.0,
    ) -> "Viewer":
        """
        Add coils to the viewer.

        Args:
            source: Coils file path or CoilsFile object
            cross_section: (width, height) for rectangular cross-section mesh
            color: Coil mesh color
            show_centerline: Also show centerline as line
            centerline_color: Color for centerline
            centerline_width: Width for centerline

        Returns:
            self for chaining
        """
        import pyvista as pv

        from .coils import CoilsFile, coils_to_mesh

        if isinstance(source, (str, Path)):
            coils = CoilsFile.from_file(source)
        else:
            coils = source

        if cross_section is not None:
            width, height = cross_section
            mesh = coils_to_mesh(coils, width=width, height=height)
            pv_mesh = pv.wrap(mesh)
            self._plotter.add_mesh(pv_mesh, color=color)
            self._meshes.append(pv_mesh)

        if show_centerline or cross_section is None:
            for fil in coils.filaments:
                pts = fil.coords
                if pts.shape[0] < 2:
                    continue
                line = pv.lines_from_points(pts)
                self._plotter.add_mesh(
                    line,
                    color=centerline_color,
                    line_width=centerline_width,
                )

        return self

    def add_vmec_lcfs(
        self,
        wout_path: Path | str,
        *,
        color: str = "red",
        opacity: float = 0.5,
        n_theta: int = 64,
        n_zeta: int = 48,
    ) -> "Viewer":
        """
        Add VMEC LCFS surface to the viewer.

        Args:
            wout_path: Path to VMEC wout file
            color: Surface color
            opacity: Surface opacity
            n_theta: Poloidal resolution
            n_zeta: Toroidal resolution

        Returns:
            self for chaining
        """
        import pyvista as pv

        from .wall import wall_mesh_from_wout

        mesh = wall_mesh_from_wout(
            wout_path, n_theta=n_theta, n_zeta=n_zeta
        )
        pv_mesh = pv.wrap(mesh)
        self._plotter.add_mesh(
            pv_mesh,
            color=color,
            opacity=opacity,
            show_edges=False,
        )
        self._meshes.append(pv_mesh)
        return self

    def add_axes(self, *, interactive: bool = True) -> "Viewer":
        """Add coordinate axes widget."""
        self._plotter.add_axes(interactive=interactive)
        return self

    def add_scalar_bar(self, title: str = "") -> "Viewer":
        """Add scalar bar for last mesh."""
        self._plotter.add_scalar_bar(title=title)
        return self

    def set_camera(
        self,
        position: tuple[float, float, float] | None = None,
        focal_point: tuple[float, float, float] | None = None,
        view_up: tuple[float, float, float] | None = None,
    ) -> "Viewer":
        """Set camera position and orientation."""
        if position is not None:
            self._plotter.camera_position = position
        if focal_point is not None:
            self._plotter.camera.focal_point = focal_point
        if view_up is not None:
            self._plotter.camera.up = view_up
        return self

    def show(self, *, interactive: bool = True) -> None:
        """
        Display the viewer.

        Args:
            interactive: Whether to allow interaction
        """
        self._plotter.show(interactive_update=interactive)

    def screenshot(
        self,
        filename: Path | str,
        *,
        window_size: tuple[int, int] = (1920, 1080),
        transparent_background: bool = False,
    ) -> Path:
        """
        Save screenshot to file.

        Args:
            filename: Output file path
            window_size: Image dimensions
            transparent_background: Use transparent background

        Returns:
            Path to saved file
        """
        filename = Path(filename)
        self._plotter.window_size = window_size
        self._plotter.screenshot(
            str(filename),
            transparent_background=transparent_background,
        )
        return filename

    def close(self) -> None:
        """Close the viewer."""
        self._plotter.close()


def quick_view(
    *sources: Path | str,
    coil_cross_section: tuple[float, float] | None = (0.4, 0.8),
    show: bool = True,
) -> Viewer:
    """
    Quick visualization of multiple geometry files.

    Automatically detects file types (STL, coils, wout) and adds them.

    Args:
        *sources: Paths to geometry files
        coil_cross_section: Cross-section for coil visualization
        show: Whether to show immediately

    Returns:
        Viewer object
    """
    viewer = Viewer()
    viewer.add_axes()

    for source in sources:
        path = Path(source)
        suffix = path.suffix.lower()
        name = path.name.lower()

        if suffix == ".stl":
            viewer.add_wall(path)
        elif suffix == ".nc" and "wout" in name:
            viewer.add_vmec_lcfs(path)
        elif "coil" in name or suffix == ".scaled":
            viewer.add_coils(path, cross_section=coil_cross_section)
        else:
            # Try to detect by content
            try:
                viewer.add_coils(path, cross_section=coil_cross_section)
            except Exception:
                try:
                    viewer.add_wall(path)
                except Exception:
                    print(f"Could not load: {path}")

    if show:
        viewer.show()

    return viewer
