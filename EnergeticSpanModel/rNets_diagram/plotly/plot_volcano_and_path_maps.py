#!/usr/bin/env python3
"""Plot volcano and minimum-path maps from kinetic_001_011_grid.csv."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np


HERE = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--grid", type=Path, default=HERE / "kinetic_001_011_grid.csv"
    )
    parser.add_argument(
        "--volcano-output", type=Path, default=HERE / "volcano_map_001_011.png"
    )
    parser.add_argument(
        "--path-output", type=Path, default=HERE / "minimum_path_map_001_011.png"
    )
    return parser.parse_args()


def load_grid(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path}: no data")

    mu_values = np.array(
        sorted({float(row["mu_MnOH2_eV"]) for row in rows}), dtype=float
    )
    water_values = np.array(
        sorted({float(row["mu_H2O_eV"]) for row in rows}), dtype=float
    )
    if len(mu_values) != len(water_values) or not np.allclose(mu_values, water_values):
        raise ValueError("This plotter expects identical square chemical-potential axes")
    expected = len(mu_values) ** 2
    if len(rows) != expected:
        raise ValueError(f"Expected {expected} grid rows, found {len(rows)}")

    row_lookup = {
        (float(row["mu_MnOH2_eV"]), float(row["mu_H2O_eV"])): row for row in rows
    }
    ordered = [
        row_lookup[(mu_mn, mu_h2o)]
        for mu_h2o in mu_values
        for mu_mn in mu_values
    ]
    shape = (len(mu_values), len(mu_values))
    data: dict[str, np.ndarray] = {}
    for surface in ("001", "011"):
        data[f"span_{surface}"] = np.array(
            [float(row[f"span_{surface}_eV"]) for row in ordered]
        ).reshape(shape)
        data[f"region_{surface}"] = np.array(
            [int(row[f"path_region_{surface}"][1:]) for row in ordered], dtype=int
        ).reshape(shape)
        data[f"cooptimal_{surface}"] = np.array(
            [int(row[f"cooptimal_paths_{surface}"]) for row in ordered], dtype=int
        ).reshape(shape)
    return mu_values, data


def nice_maximum(value: float) -> float:
    if value <= 1:
        return math.ceil(value * 10) / 10
    return math.ceil(value * 2) / 2


def add_extreme_markers(ax: plt.Axes, values: np.ndarray, mu_values: np.ndarray) -> None:
    minimum_index = np.unravel_index(int(np.argmin(values)), values.shape)
    maximum_index = np.unravel_index(int(np.argmax(values)), values.shape)
    minimum = float(values[minimum_index])
    maximum = float(values[maximum_index])
    ax.scatter(
        mu_values[minimum_index[1]],
        mu_values[minimum_index[0]],
        marker="*",
        s=150,
        facecolor="white",
        edgecolor="black",
        linewidth=1.2,
        zorder=6,
        label=rf"min $\delta E$ = {minimum:.2f} eV",
    )
    ax.scatter(
        mu_values[maximum_index[1]],
        mu_values[maximum_index[0]],
        marker="*",
        s=150,
        facecolor="#ff4057",
        edgecolor="white",
        linewidth=0.8,
        zorder=6,
        label=rf"max $\delta E$ = {maximum:.2f} eV",
    )
    ax.legend(loc="best", fontsize=8, framealpha=0.9)


def plot_volcano(
    output: Path, mu_values: np.ndarray, data: dict[str, np.ndarray]
) -> None:
    mu_mn, mu_h2o = np.meshgrid(mu_values, mu_values, indexing="xy")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.3), constrained_layout=True)
    for ax, surface, subtitle in (
        (axes[0], "001", "2-layer repeat; 2.354 Å"),
        (axes[1], "011", "4-layer repeat; 4.906 Å"),
    ):
        values = data[f"span_{surface}"]
        maximum = nice_maximum(float(values.max()))
        filled_levels = np.linspace(0.0, maximum, 22)
        contour_levels = np.linspace(0.0, maximum, 9)[1:-1]
        filled = ax.contourf(
            mu_mn,
            mu_h2o,
            values,
            levels=filled_levels,
            cmap="viridis_r",
            extend="max",
        )
        lines = ax.contour(
            mu_mn,
            mu_h2o,
            values,
            levels=contour_levels,
            colors="#404040",
            linewidths=0.55,
            alpha=0.7,
        )
        ax.clabel(lines, inline=True, fontsize=7, fmt="%.2f")
        add_extreme_markers(ax, values, mu_values)
        ax.set_title(f"[{surface}] energetic-span volcano\n{subtitle}")
        ax.set_xlabel(r"$\Delta\mu_{\mathrm{Mn(OH)_2}}$ (eV)")
        ax.set_ylabel(r"$\Delta\mu_{\mathrm{H_2O}}$ (eV)")
        ax.set_aspect("equal")
        colorbar = fig.colorbar(filled, ax=ax, shrink=0.9)
        colorbar.set_label(r"minimum $\delta E$ (eV)")
    fig.suptitle(
        "Minimum energetic-span volcano maps",
        fontsize=15,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def boundary_segments(region_grid: np.ndarray, mu_values: np.ndarray) -> list[list[tuple[float, float]]]:
    step = float(mu_values[1] - mu_values[0])
    segments: list[list[tuple[float, float]]] = []
    rows, columns = region_grid.shape
    for row in range(rows):
        y = mu_values[row]
        for column in range(columns - 1):
            if region_grid[row, column] == region_grid[row, column + 1]:
                continue
            x = 0.5 * (mu_values[column] + mu_values[column + 1])
            segments.append([(x, y - step / 2), (x, y + step / 2)])
    for row in range(rows - 1):
        y = 0.5 * (mu_values[row] + mu_values[row + 1])
        for column in range(columns):
            if region_grid[row, column] == region_grid[row + 1, column]:
                continue
            x = mu_values[column]
            segments.append([(x - step / 2, y), (x + step / 2, y)])
    return segments


def connected_components(region_grid: np.ndarray) -> list[tuple[int, list[tuple[int, int]]]]:
    rows, columns = region_grid.shape
    visited = np.zeros_like(region_grid, dtype=bool)
    components = []
    for row in range(rows):
        for column in range(columns):
            if visited[row, column]:
                continue
            value = int(region_grid[row, column])
            stack = [(row, column)]
            visited[row, column] = True
            cells = []
            while stack:
                current = stack.pop()
                cells.append(current)
                for neighbor in (
                    (current[0] - 1, current[1]),
                    (current[0] + 1, current[1]),
                    (current[0], current[1] - 1),
                    (current[0], current[1] + 1),
                ):
                    if not (0 <= neighbor[0] < rows and 0 <= neighbor[1] < columns):
                        continue
                    if visited[neighbor] or int(region_grid[neighbor]) != value:
                        continue
                    visited[neighbor] = True
                    stack.append(neighbor)
            components.append((value, cells))
    return components


def annotate_regions(ax: plt.Axes, region_grid: np.ndarray, mu_values: np.ndarray) -> None:
    minimum_cells = max(5, int(region_grid.size * 0.002))
    for region, cells in connected_components(region_grid):
        if len(cells) < minimum_cells:
            continue
        mean_row = sum(cell[0] for cell in cells) / len(cells)
        mean_column = sum(cell[1] for cell in cells) / len(cells)
        row, column = min(
            cells,
            key=lambda cell: (cell[0] - mean_row) ** 2 + (cell[1] - mean_column) ** 2,
        )
        ax.text(
            mu_values[column],
            mu_values[row],
            f"R{region:03d}",
            ha="center",
            va="center",
            fontsize=6.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.8},
            zorder=5,
        )


def plot_path_map(
    output: Path, mu_values: np.ndarray, data: dict[str, np.ndarray]
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.3, 5.4), constrained_layout=True)
    extent = [
        mu_values[0],
        mu_values[-1],
        mu_values[0],
        mu_values[-1],
    ]
    for ax, surface, subtitle in (
        (axes[0], "001", "2-layer minimum path"),
        (axes[1], "011", "4-layer minimum path"),
    ):
        regions = data[f"region_{surface}"]
        active = np.array(sorted(np.unique(regions)), dtype=int)
        rank = {region: index for index, region in enumerate(active)}
        ranked = np.vectorize(rank.__getitem__)(regions)
        colors = plt.colormaps["turbo"](np.linspace(0.03, 0.97, len(active)))
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(np.arange(len(active) + 1) - 0.5, len(active))
        ax.imshow(
            ranked,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
        )
        segments = boundary_segments(regions, mu_values)
        if segments:
            ax.add_collection(
                LineCollection(segments, colors="black", linewidths=0.65, alpha=0.9)
            )
        annotate_regions(ax, regions, mu_values)
        cooptimal_fraction = float(np.mean(data[f"cooptimal_{surface}"] > 1))
        maximum_cooptimal = int(data[f"cooptimal_{surface}"].max())
        ax.set_title(
            f"[{surface}] {subtitle}\n"
            f"co-optimal grid: {100 * cooptimal_fraction:.1f}% "
            f"(maximum {maximum_cooptimal} paths)"
        )
        ax.set_xlabel(r"$\Delta\mu_{\mathrm{Mn(OH)_2}}$ (eV)")
        ax.set_ylabel(r"$\Delta\mu_{\mathrm{H_2O}}$ (eV)")
        ax.set_xlim(mu_values[0], mu_values[-1])
        ax.set_ylim(mu_values[0], mu_values[-1])
        ax.set_aspect("equal")
    fig.suptitle(
        "Representative minimum energetic-span path regions\n"
        "Black lines = path-change boundaries; R-codes map to kinetic_001_011_path_regions.csv",
        fontsize=14,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    mu_values, data = load_grid(args.grid)
    plot_volcano(args.volcano_output, mu_values, data)
    plot_path_map(args.path_output, mu_values, data)
    print(f"Wrote {args.volcano_output}")
    print(f"Wrote {args.path_output}")


if __name__ == "__main__":
    main()
