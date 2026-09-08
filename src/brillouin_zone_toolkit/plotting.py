"""Matplotlib visualizations for zones and symmetry-reduced point clouds."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from .brillouin import BrillouinZone
from .lattices import Lattice
from .symmetry import lattice_symmetry_operations, reduce_points


def _axis_limits(zone: BrillouinZone, padding: float = 0.15) -> tuple[np.ndarray, np.ndarray]:
    low = np.min(zone.vertices, axis=0)
    high = np.max(zone.vertices, axis=0)
    span = np.maximum(high - low, 1e-12)
    return low - padding * span, high + padding * span


def _set_equal_3d_axes(ax: Any, zone: BrillouinZone) -> None:
    low, high = _axis_limits(zone)
    ax.set_xlim(float(low[0]), float(high[0]))
    ax.set_ylim(float(low[1]), float(high[1]))
    ax.set_zlim(float(low[2]), float(high[2]))
    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect(tuple(high - low))
    ax.set_xlabel(r"$k_x$")
    ax.set_ylabel(r"$k_y$")
    ax.set_zlabel(r"$k_z$")


def plot_brillouin_zone(
    zone: BrillouinZone,
    *,
    ax: Any = None,
    face_color: str = "#4C78A8",
    edge_color: str = "#17324D",
    alpha: float = 0.35,
    show_vertices: bool = True,
    title: Optional[str] = None,
) -> Any:
    """Plot a zone as a translucent 3-D convex polyhedron.

    Matplotlib is imported lazily so reciprocal-lattice calculations remain
    usable in headless or minimal Python environments.
    """

    if not isinstance(zone, BrillouinZone):
        raise TypeError("zone must be a BrillouinZone")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    polygons = [zone.vertices[list(face)] for face in zone.faces]
    collection = Poly3DCollection(
        polygons, facecolors=face_color, edgecolors=edge_color, linewidths=0.8, alpha=alpha
    )
    ax.add_collection3d(collection)
    if show_vertices:
        ax.scatter(
            zone.vertices[:, 0],
            zone.vertices[:, 1],
            zone.vertices[:, 2],
            color=edge_color,
            s=12,
            depthshade=False,
        )
    ax.scatter([0.0], [0.0], [0.0], color="#D62728", s=28, marker="+")
    _set_equal_3d_axes(ax, zone)
    if title:
        ax.set_title(title)
    return ax


def plot_symmetry_reduction(
    zone: BrillouinZone,
    points: Any,
    *,
    lattice: Optional[Lattice] = None,
    operations: Any = None,
    ax: Any = None,
    all_color: str = "#B8C6D9",
    reduced_color: str = "#F58518",
    title: Optional[str] = None,
) -> Any:
    """Plot an input point cloud and one representative from each orbit."""

    if operations is None:
        if lattice is None:
            raise ValueError("provide lattice or operations")
        operations = lattice_symmetry_operations(lattice)
    points_array = np.asarray(points, dtype=float)
    if points_array.ndim != 2 or points_array.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    representatives = reduce_points(points_array, operations)
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    plot_brillouin_zone(zone, ax=ax, alpha=0.10, show_vertices=False)
    ax.scatter(
        points_array[:, 0], points_array[:, 1], points_array[:, 2],
        color=all_color, s=5, depthshade=False, label="sampled points",
    )
    ax.scatter(
        representatives[:, 0], representatives[:, 1], representatives[:, 2],
        color=reduced_color, s=14, depthshade=False, label="symmetry representatives",
    )
    ax.legend(loc="upper left", fontsize="small")
    if title:
        ax.set_title(title)
    return ax


def sample_cartesian_grid(zone: BrillouinZone, points_per_axis: int = 15) -> np.ndarray:
    """Sample a deterministic Cartesian grid and retain points in ``zone``."""

    if not isinstance(points_per_axis, (int, np.integer)) or points_per_axis < 2:
        raise ValueError("points_per_axis must be an integer >= 2")
    low, high = _axis_limits(zone, padding=0.0)
    axes = [np.linspace(low[index], high[index], int(points_per_axis)) for index in range(3)]
    mesh = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape((-1, 3))
    return mesh[zone.contains(mesh)]


def plot_reciprocal_lattice(
    reciprocal_points: Any, *, ax: Any = None, title: Optional[str] = None
) -> Any:
    """Plot reciprocal-lattice points as a 3-D scatter plot."""

    points = np.asarray(reciprocal_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("reciprocal_points must have shape (n, 3)")
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=9, alpha=0.7)
    ax.scatter([0.0], [0.0], [0.0], color="#D62728", s=28, marker="+")
    ax.set_xlabel(r"$G_x$")
    ax.set_ylabel(r"$G_y$")
    ax.set_zlabel(r"$G_z$")
    if title:
        ax.set_title(title)
    return ax
