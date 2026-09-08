#!/usr/bin/env python3
"""Reproduce the figures shipped with the Brillouin-zone toolkit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from brillouin_zone_toolkit import (
    bcc,
    fcc,
    first_brillouin_zone,
    hexagonal,
    simple_cubic,
)


def _matplotlib():
    """Import Matplotlib lazily and select a headless backend."""

    try:
        import matplotlib
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise RuntimeError(
            "figure reproduction requires Matplotlib; install the 'plot' extra"
        ) from error
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt
from brillouin_zone_toolkit.plotting import (
    plot_brillouin_zone,
    plot_symmetry_reduction,
    sample_cartesian_grid,
)
from brillouin_zone_toolkit.symmetry import lattice_symmetry_operations, reduce_points


def _save(fig: object, path: Path) -> None:
    plt = _matplotlib()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def generate_figures(output: Path) -> dict[str, dict[str, float | int | str]]:
    plt = _matplotlib()
    output.mkdir(parents=True, exist_ok=True)
    lattices = [simple_cubic(), fcc(), bcc(), hexagonal()]
    zones = [(lattice, first_brillouin_zone(lattice, search=2)) for lattice in lattices]
    summary: dict[str, dict[str, float | int | str]] = {}

    for lattice, zone in zones:
        slug = lattice.name.split(" ")[0].lower()
        if slug == "face-centred":
            slug = "fcc"
        elif slug == "body-centred":
            slug = "bcc"
        elif slug == "simple":
            slug = "sc"
        figure = plt.figure(figsize=(6.4, 5.2))
        axis = figure.add_subplot(111, projection="3d")
        plot_brillouin_zone(zone, ax=axis, title=f"First Brillouin zone — {lattice.name}")
        _save(figure, output / f"{slug}-first-brillouin-zone.png")
        summary[slug] = {
            "lattice": lattice.name,
            "vertices": zone.n_vertices,
            "faces": zone.n_faces,
            "volume": zone.volume,
            "reciprocal_cell_volume": lattice.reciprocal_volume,
        }

    figure = plt.figure(figsize=(11, 8))
    for index, (lattice, zone) in enumerate(zones, start=1):
        axis = figure.add_subplot(2, 2, index, projection="3d")
        plot_brillouin_zone(zone, ax=axis, title=lattice.name, show_vertices=False)
    figure.suptitle("Gallery: first Brillouin zones")
    _save(figure, output / "gallery.png")

    fcc_lattice = fcc()
    fcc_zone = first_brillouin_zone(fcc_lattice, search=2)
    points = sample_cartesian_grid(fcc_zone, points_per_axis=17)
    operations = lattice_symmetry_operations(fcc_lattice)
    reduced = reduce_points(points, operations)
    figure = plt.figure(figsize=(7.2, 5.8))
    axis = figure.add_subplot(111, projection="3d")
    plot_symmetry_reduction(
        fcc_zone,
        points,
        operations=operations,
        ax=axis,
        title=f"FCC symmetry reduction ({len(points)} → {len(reduced)} representatives)",
    )
    _save(figure, output / "fcc-symmetry-reduction.png")
    summary["fcc-symmetry-reduction"] = {
        "sampled_points": len(points),
        "representatives": len(reduced),
        "operations": len(operations),
    }

    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs" / "figures",
        help="directory receiving PNG figures and summary.json",
    )
    arguments = parser.parse_args()
    summary = generate_figures(arguments.output)
    for name, values in summary.items():
        print(f"{name}: {values}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
