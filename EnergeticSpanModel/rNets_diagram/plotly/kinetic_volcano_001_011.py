#!/usr/bin/env python3
"""Compare [001] and [011] growth kinetics on a chemical-potential grid.

The calculation follows the minimum energetic-span convention already used by
the Plotly notebooks in this directory.  One repeat cycle is S000 -> S322 for
[001] and S000 -> S644 for [011].  Cycle frequencies are converted to normal
growth velocities using repeat heights of 2.354 and 4.906 Angstrom,
respectively.

The energy CSV files currently contain relaxed states rather than explicit
transition states.  Consequently, the reported TOFs and velocities are
energetic-span-derived kinetic proxies, not a full microkinetic model.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
import numpy as np


HERE = Path(__file__).resolve().parent
KB_EV_PER_K = 8.617333262145e-5
KB_J_PER_K = 1.380649e-23
PLANCK_J_S = 6.62607015e-34
TIE_TOLERANCE_EV = 1e-6


@dataclass(frozen=True)
class SurfaceSpec:
    label: str
    csv_file: Path
    start_state: str
    final_state: str
    repeat_height_angstrom: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--energy-001", type=Path, default=HERE / "energy_001_allmag.csv"
    )
    parser.add_argument(
        "--energy-011", type=Path, default=HERE / "energy_011_1234.csv"
    )
    parser.add_argument("--height-001", type=float, default=2.354)
    parser.add_argument("--height-011", type=float, default=4.906)
    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--mu-min", type=float, default=-4.0)
    parser.add_argument("--mu-max", type=float, default=0.0)
    parser.add_argument("--mu-step", type=float, default=0.05)
    parser.add_argument("--path-batch", type=int, default=4000)
    parser.add_argument("--grid-batch", type=int, default=243)
    parser.add_argument(
        "--output-prefix", type=Path, default=HERE / "kinetic_001_011"
    )
    return parser.parse_args()


def load_states(path: Path) -> list[dict]:
    states = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"state", "coordinate", "relative_E"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row in reader:
            name = row["state"].strip()
            energy = row["relative_E"].strip()
            if not name or not energy:
                continue
            if len(name) != 4 or not name.startswith("S") or not name[1:].isdigit():
                raise ValueError(f"{path}: invalid state label {name!r}")
            o_count, c_count, x_count = map(int, name[1:])
            states.append(
                {
                    "name": name,
                    "coordinate": float(row["coordinate"]),
                    "E0": float(energy),
                    "O": o_count,
                    "C": c_count,
                    "X": x_count,
                }
            )
    states.sort(key=lambda state: (state["coordinate"], state["name"]))
    names = [state["name"] for state in states]
    if len(names) != len(set(names)):
        raise ValueError(f"{path}: duplicate state labels are not supported")
    return states


def classify_transition(initial: dict, final: dict) -> str | None:
    delta = (
        final["O"] - initial["O"],
        final["C"] - initial["C"],
        final["X"] - initial["X"],
    )
    return {
        (1, 0, 0): "OL",
        (0, 1, 0): "C",
        (0, 0, 1): "SR",
        (1, 1, 0): "OX",
    }.get(delta)


def build_edges(states: list[dict]) -> list[tuple[str, str, str]]:
    edges = []
    for initial in states:
        for final in states:
            if final["coordinate"] <= initial["coordinate"]:
                continue
            mechanism = classify_transition(initial, final)
            if mechanism is None:
                continue
            expected_step = 2 if mechanism == "OX" else 1
            if math.isclose(
                final["coordinate"] - initial["coordinate"],
                expected_step,
                abs_tol=1e-12,
            ):
                edges.append((initial["name"], final["name"], mechanism))
    return edges


def enumerate_paths(
    edges: list[tuple[str, str, str]], start_state: str, final_state: str
) -> list[dict]:
    adjacency: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for start, end, mechanism in edges:
        adjacency[start].append((end, mechanism))
    for name in adjacency:
        adjacency[name].sort()

    paths = []

    def dfs(node: str, states: list[str], mechanisms: list[str]) -> None:
        if node == final_state:
            paths.append({"states": states.copy(), "mechanisms": mechanisms.copy()})
            return
        for next_node, mechanism in adjacency.get(node, []):
            if next_node not in states:
                dfs(next_node, states + [next_node], mechanisms + [mechanism])

    dfs(start_state, [start_state], [])
    paths.sort(key=lambda path: " -> ".join(path["states"]))
    for index, path in enumerate(paths, start=1):
        path["path_id"] = f"P{index:05d}"
    if not paths:
        raise ValueError(f"No complete path from {start_state} to {final_state}")
    return paths


def path_span_matrix(
    candidate_indices: np.ndarray,
    candidate_valid: np.ndarray,
    state_energies: np.ndarray,
    grid_indices: np.ndarray,
    reaction_energy: np.ndarray,
) -> np.ndarray:
    """Return energetic spans for path batch x chemical-potential batch."""
    batch, length = candidate_indices.shape
    grids = len(grid_indices)
    minimum_prefix = np.full((batch, grids), np.inf)
    maximum_prefix = np.full((batch, grids), -np.inf)
    maximum_rise = np.zeros((batch, grids))
    maximum_drop = np.zeros((batch, grids))

    for column in range(length):
        valid = candidate_valid[:, column][:, None]
        indices = candidate_indices[:, column]
        safe_indices = np.where(indices >= 0, indices, 0)
        energy = state_energies[safe_indices[:, None], grid_indices[None, :]]
        new_minimum = np.minimum(minimum_prefix, energy)
        maximum_rise = np.where(
            valid, np.maximum(maximum_rise, energy - new_minimum), maximum_rise
        )
        maximum_drop = np.where(
            valid, np.maximum(maximum_drop, maximum_prefix - energy), maximum_drop
        )
        minimum_prefix = np.where(valid, new_minimum, minimum_prefix)
        maximum_prefix = np.where(
            valid, np.maximum(maximum_prefix, energy), maximum_prefix
        )

    return np.maximum(
        maximum_rise,
        maximum_drop + reaction_energy[grid_indices][None, :],
    )


def prepare_surface(spec: SurfaceSpec) -> dict:
    states = load_states(spec.csv_file)
    state_index = {state["name"]: index for index, state in enumerate(states)}
    for endpoint in (spec.start_state, spec.final_state):
        if endpoint not in state_index:
            raise ValueError(f"{spec.csv_file}: endpoint {endpoint} missing")

    edges = build_edges(states)
    paths = enumerate_paths(edges, spec.start_state, spec.final_state)
    maximum_candidates = max(len(path["states"]) - 1 for path in paths)
    candidate_indices = np.full(
        (len(paths), maximum_candidates), -1, dtype=np.int16
    )
    for row, path in enumerate(paths):
        names = path["states"][:-1]
        candidate_indices[row, : len(names)] = [state_index[name] for name in names]

    coefficients = np.array(
        [
            [state["E0"], -state["O"], state["C"] + 0.5 * state["X"]]
            for state in states
        ],
        dtype=float,
    )
    print(
        f"{spec.label}: {len(states)} states, {len(edges)} edges, "
        f"{len(paths)} complete paths",
        flush=True,
    )
    return {
        "spec": spec,
        "states": states,
        "state_index": state_index,
        "paths": paths,
        "candidate_indices": candidate_indices,
        "candidate_valid": candidate_indices >= 0,
        "coefficients": coefficients,
    }


def evaluate_surface(
    prepared: dict,
    mu_mn_flat: np.ndarray,
    mu_h2o_flat: np.ndarray,
    path_batch: int,
    grid_batch: int,
) -> dict:
    condition = np.vstack(
        [np.ones_like(mu_mn_flat), mu_mn_flat, mu_h2o_flat]
    )
    state_energies = prepared["coefficients"] @ condition
    state_index = prepared["state_index"]
    spec = prepared["spec"]
    reaction_energy = (
        state_energies[state_index[spec.final_state]]
        - state_energies[state_index[spec.start_state]]
    )

    grid_count = len(mu_mn_flat)
    best_span = np.full(grid_count, np.inf)
    best_path_index = np.full(grid_count, -1, dtype=np.int32)
    path_count = len(prepared["paths"])

    for grid_start in range(0, grid_count, grid_batch):
        grid_indices = np.arange(grid_start, min(grid_start + grid_batch, grid_count))
        local_span = np.full(len(grid_indices), np.inf)
        local_path = np.full(len(grid_indices), -1, dtype=np.int32)
        for path_start in range(0, path_count, path_batch):
            path_end = min(path_start + path_batch, path_count)
            spans = path_span_matrix(
                prepared["candidate_indices"][path_start:path_end],
                prepared["candidate_valid"][path_start:path_end],
                state_energies,
                grid_indices,
                reaction_energy,
            )
            indices = np.argmin(spans, axis=0)
            values = spans[indices, np.arange(len(grid_indices))]
            better = values < local_span - TIE_TOLERANCE_EV
            local_span[better] = values[better]
            local_path[better] = path_start + indices[better]
        best_span[grid_indices] = local_span
        best_path_index[grid_indices] = local_path
        print(
            f"{spec.label}: evaluated {grid_indices[-1] + 1}/{grid_count} grid points",
            flush=True,
        )

    cooptimal_path_count = np.zeros(grid_count, dtype=np.int32)
    for grid_start in range(0, grid_count, grid_batch):
        grid_indices = np.arange(grid_start, min(grid_start + grid_batch, grid_count))
        local_count = np.zeros(len(grid_indices), dtype=np.int32)
        for path_start in range(0, path_count, path_batch):
            path_end = min(path_start + path_batch, path_count)
            spans = path_span_matrix(
                prepared["candidate_indices"][path_start:path_end],
                prepared["candidate_valid"][path_start:path_end],
                state_energies,
                grid_indices,
                reaction_energy,
            )
            local_count += np.sum(
                np.abs(spans - best_span[grid_indices][None, :])
                <= TIE_TOLERANCE_EV,
                axis=0,
            ).astype(np.int32)
        cooptimal_path_count[grid_indices] = local_count

    tdi = []
    tdts = []
    for grid, path_index in enumerate(best_path_index):
        names = prepared["paths"][path_index]["states"][:-1]
        indices = np.array([state_index[name] for name in names])
        energies = state_energies[indices, grid]
        span_matrix = energies[:, None] - energies[None, :]
        row_index = np.arange(len(names))[:, None]
        column_index = np.arange(len(names))[None, :]
        span_matrix = span_matrix + np.where(
            row_index < column_index, reaction_energy[grid], 0.0
        )
        tdts_index, tdi_index = np.unravel_index(
            int(np.argmax(span_matrix)), span_matrix.shape
        )
        tdts.append(names[tdts_index])
        tdi.append(names[tdi_index])

    return {
        "best_span": best_span,
        "best_path_index": best_path_index,
        "reaction_energy": reaction_energy,
        "tdi": np.array(tdi, dtype=object),
        "tdts": np.array(tdts, dtype=object),
        "cooptimal_path_count": cooptimal_path_count,
        "state_energies": state_energies,
    }


def assign_region_labels(prepared: dict, result: dict) -> dict[int, str]:
    active = sorted(set(map(int, result["best_path_index"])))
    return {path_index: f"R{rank:03d}" for rank, path_index in enumerate(active, 1)}


def boundary_segments(
    region_grid: np.ndarray, mu_values: np.ndarray
) -> list[tuple[tuple[float, float], tuple[float, float], int, int]]:
    step = float(mu_values[1] - mu_values[0]) if len(mu_values) > 1 else 0.0
    segments = []
    rows, columns = region_grid.shape

    for row in range(rows):
        y = mu_values[row]
        for column in range(columns - 1):
            left = int(region_grid[row, column])
            right = int(region_grid[row, column + 1])
            if left == right:
                continue
            x = 0.5 * (mu_values[column] + mu_values[column + 1])
            segments.append(((x, y - step / 2), (x, y + step / 2), left, right))

    for row in range(rows - 1):
        y = 0.5 * (mu_values[row] + mu_values[row + 1])
        for column in range(columns):
            lower = int(region_grid[row, column])
            upper = int(region_grid[row + 1, column])
            if lower == upper:
                continue
            x = mu_values[column]
            segments.append(((x - step / 2, y), (x + step / 2, y), lower, upper))
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
            queue = deque([(row, column)])
            visited[row, column] = True
            cells = []
            while queue:
                current_row, current_column = queue.popleft()
                cells.append((current_row, current_column))
                for next_row, next_column in (
                    (current_row - 1, current_column),
                    (current_row + 1, current_column),
                    (current_row, current_column - 1),
                    (current_row, current_column + 1),
                ):
                    if not (0 <= next_row < rows and 0 <= next_column < columns):
                        continue
                    if visited[next_row, next_column]:
                        continue
                    if int(region_grid[next_row, next_column]) != value:
                        continue
                    visited[next_row, next_column] = True
                    queue.append((next_row, next_column))
            components.append((value, cells))
    return components


def draw_path_boundaries_and_labels(
    ax: plt.Axes,
    region_grid: np.ndarray,
    mu_values: np.ndarray,
    region_labels: dict[int, str],
) -> None:
    segments = boundary_segments(region_grid, mu_values)
    if segments:
        ax.add_collection(
            LineCollection(
                [[start, end] for start, end, _, _ in segments],
                colors="black",
                linewidths=0.8,
                alpha=0.9,
            )
        )

    minimum_label_cells = max(4, int(region_grid.size * 0.002))
    for region, cells in connected_components(region_grid):
        if len(cells) < minimum_label_cells:
            continue
        mean_row = int(round(sum(cell[0] for cell in cells) / len(cells)))
        mean_column = int(round(sum(cell[1] for cell in cells) / len(cells)))
        candidates = sorted(
            cells,
            key=lambda cell: (cell[0] - mean_row) ** 2 + (cell[1] - mean_column) ** 2,
        )
        center_row, center_column = candidates[0]
        ax.text(
            mu_values[center_column],
            mu_values[center_row],
            region_labels[region],
            ha="center",
            va="center",
            fontsize=7,
            color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 1},
        )


def write_grid_csv(
    path: Path,
    mu_mn: np.ndarray,
    mu_h2o: np.ndarray,
    prepared_001: dict,
    result_001: dict,
    labels_001: dict[int, str],
    prepared_011: dict,
    result_011: dict,
    labels_011: dict[int, str],
    temperature: float,
) -> None:
    thermal_ev = KB_EV_PER_K * temperature
    log_prefactor = math.log10(KB_J_PER_K * temperature / PLANCK_J_S)
    log_tof_001 = log_prefactor - result_001["best_span"] / (thermal_ev * math.log(10))
    log_tof_011 = log_prefactor - result_011["best_span"] / (thermal_ev * math.log(10))
    height_001 = prepared_001["spec"].repeat_height_angstrom
    height_011 = prepared_011["spec"].repeat_height_angstrom
    log_velocity_001 = math.log10(height_001) + log_tof_001
    log_velocity_011 = math.log10(height_011) + log_tof_011

    columns = [
        "mu_MnOH2_eV",
        "mu_H2O_eV",
        "span_001_eV",
        "span_011_eV",
        "delta_span_011_minus_001_eV",
        "log10_TOF_001_s-1",
        "log10_TOF_011_s-1",
        "log10_velocity_001_A_s-1",
        "log10_velocity_011_A_s-1",
        "log10_velocity_ratio_001_over_011",
        "path_region_001",
        "path_id_001",
        "TDI_001",
        "TDTS_001",
        "cooptimal_paths_001",
        "path_region_011",
        "path_id_011",
        "TDI_011",
        "TDTS_011",
        "cooptimal_paths_011",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for grid in range(len(mu_mn)):
            index_001 = int(result_001["best_path_index"][grid])
            index_011 = int(result_011["best_path_index"][grid])
            writer.writerow(
                {
                    "mu_MnOH2_eV": f"{mu_mn[grid]:.8g}",
                    "mu_H2O_eV": f"{mu_h2o[grid]:.8g}",
                    "span_001_eV": f"{result_001['best_span'][grid]:.9g}",
                    "span_011_eV": f"{result_011['best_span'][grid]:.9g}",
                    "delta_span_011_minus_001_eV": f"{result_011['best_span'][grid] - result_001['best_span'][grid]:.9g}",
                    "log10_TOF_001_s-1": f"{log_tof_001[grid]:.9g}",
                    "log10_TOF_011_s-1": f"{log_tof_011[grid]:.9g}",
                    "log10_velocity_001_A_s-1": f"{log_velocity_001[grid]:.9g}",
                    "log10_velocity_011_A_s-1": f"{log_velocity_011[grid]:.9g}",
                    "log10_velocity_ratio_001_over_011": f"{log_velocity_001[grid] - log_velocity_011[grid]:.9g}",
                    "path_region_001": labels_001[index_001],
                    "path_id_001": prepared_001["paths"][index_001]["path_id"],
                    "TDI_001": result_001["tdi"][grid],
                    "TDTS_001": result_001["tdts"][grid],
                    "cooptimal_paths_001": int(result_001["cooptimal_path_count"][grid]),
                    "path_region_011": labels_011[index_011],
                    "path_id_011": prepared_011["paths"][index_011]["path_id"],
                    "TDI_011": result_011["tdi"][grid],
                    "TDTS_011": result_011["tdts"][grid],
                    "cooptimal_paths_011": int(result_011["cooptimal_path_count"][grid]),
                }
            )


def write_region_csv(
    path: Path,
    surfaces: list[tuple[dict, dict, dict[int, str]]],
    mu_mn: np.ndarray,
    mu_h2o: np.ndarray,
) -> None:
    columns = [
        "surface",
        "region",
        "path_id",
        "grid_points",
        "mu_MnOH2_min_eV",
        "mu_MnOH2_max_eV",
        "mu_H2O_min_eV",
        "mu_H2O_max_eV",
        "states",
        "mechanisms",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for prepared, result, labels in surfaces:
            for path_index, region in labels.items():
                selected = result["best_path_index"] == path_index
                path_info = prepared["paths"][path_index]
                writer.writerow(
                    {
                        "surface": prepared["spec"].label,
                        "region": region,
                        "path_id": path_info["path_id"],
                        "grid_points": int(selected.sum()),
                        "mu_MnOH2_min_eV": f"{mu_mn[selected].min():.8g}",
                        "mu_MnOH2_max_eV": f"{mu_mn[selected].max():.8g}",
                        "mu_H2O_min_eV": f"{mu_h2o[selected].min():.8g}",
                        "mu_H2O_max_eV": f"{mu_h2o[selected].max():.8g}",
                        "states": " -> ".join(path_info["states"]),
                        "mechanisms": " -> ".join(path_info["mechanisms"]),
                    }
                )


def write_boundary_csv(
    path: Path,
    surfaces: list[tuple[dict, dict, dict[int, str]]],
    mu_values: np.ndarray,
    grid_shape: tuple[int, int],
) -> None:
    columns = [
        "surface",
        "x1_mu_MnOH2_eV",
        "y1_mu_H2O_eV",
        "x2_mu_MnOH2_eV",
        "y2_mu_H2O_eV",
        "region_a",
        "region_b",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for prepared, result, labels in surfaces:
            region_grid = result["best_path_index"].reshape(grid_shape)
            for start, end, first, second in boundary_segments(region_grid, mu_values):
                writer.writerow(
                    {
                        "surface": prepared["spec"].label,
                        "x1_mu_MnOH2_eV": f"{start[0]:.8g}",
                        "y1_mu_H2O_eV": f"{start[1]:.8g}",
                        "x2_mu_MnOH2_eV": f"{end[0]:.8g}",
                        "y2_mu_H2O_eV": f"{end[1]:.8g}",
                        "region_a": labels[first],
                        "region_b": labels[second],
                    }
                )


def linear_zero(first_coordinate: float, second_coordinate: float, first: float, second: float) -> float:
    if abs(second - first) < 1e-15:
        return 0.5 * (first_coordinate + second_coordinate)
    fraction = -first / (second - first)
    return first_coordinate + fraction * (second_coordinate - first_coordinate)


def write_equal_growth_csv(
    path: Path,
    mu_values: np.ndarray,
    log_velocity_ratio: np.ndarray,
) -> None:
    """Write linearly interpolated zero crossings along both grid directions."""
    columns = [
        "scan_axis",
        "fixed_chemical_potential_eV",
        "equal_growth_crossing_eV",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for column, fixed_mu_mn in enumerate(mu_values):
            values = log_velocity_ratio[:, column]
            for row in range(len(mu_values) - 1):
                if values[row] * values[row + 1] > 0:
                    continue
                crossing = linear_zero(
                    mu_values[row], mu_values[row + 1], values[row], values[row + 1]
                )
                writer.writerow(
                    {
                        "scan_axis": "mu_H2O_at_fixed_mu_MnOH2",
                        "fixed_chemical_potential_eV": f"{fixed_mu_mn:.8g}",
                        "equal_growth_crossing_eV": f"{crossing:.9g}",
                    }
                )
        for row, fixed_mu_h2o in enumerate(mu_values):
            values = log_velocity_ratio[row, :]
            for column in range(len(mu_values) - 1):
                if values[column] * values[column + 1] > 0:
                    continue
                crossing = linear_zero(
                    mu_values[column],
                    mu_values[column + 1],
                    values[column],
                    values[column + 1],
                )
                writer.writerow(
                    {
                        "scan_axis": "mu_MnOH2_at_fixed_mu_H2O",
                        "fixed_chemical_potential_eV": f"{fixed_mu_h2o:.8g}",
                        "equal_growth_crossing_eV": f"{crossing:.9g}",
                    }
                )


def make_figure(
    path: Path,
    mu_values: np.ndarray,
    prepared_001: dict,
    result_001: dict,
    labels_001: dict[int, str],
    prepared_011: dict,
    result_011: dict,
    labels_011: dict[int, str],
    temperature: float,
) -> None:
    shape = (len(mu_values), len(mu_values))
    span_001 = result_001["best_span"].reshape(shape)
    span_011 = result_011["best_span"].reshape(shape)
    region_001 = result_001["best_path_index"].reshape(shape)
    region_011 = result_011["best_path_index"].reshape(shape)

    thermal_ev = KB_EV_PER_K * temperature
    log_velocity_ratio = (
        math.log10(
            prepared_001["spec"].repeat_height_angstrom
            / prepared_011["spec"].repeat_height_angstrom
        )
        + (span_011 - span_001) / (thermal_ev * math.log(10))
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    extent = [mu_values[0], mu_values[-1], mu_values[0], mu_values[-1]]
    shared_min = float(min(span_001.min(), span_011.min()))
    shared_max = float(max(span_001.max(), span_011.max()))

    for ax, span, regions, labels, title in (
        (axes[0], span_001, region_001, labels_001, "[001]: 2-layer repeat (2.354 Å)"),
        (axes[1], span_011, region_011, labels_011, "[011]: 4-layer repeat (4.906 Å)"),
    ):
        image = ax.imshow(
            span,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="viridis",
            vmin=shared_min,
            vmax=shared_max,
            interpolation="nearest",
        )
        draw_path_boundaries_and_labels(ax, regions, mu_values, labels)
        ax.set_title(title)
        ax.set_xlabel(r"$\Delta\mu_{\mathrm{Mn(OH)_2}}$ (eV)")
        ax.set_ylabel(r"$\Delta\mu_{\mathrm{H_2O}}$ (eV)")
        colorbar = fig.colorbar(image, ax=ax, shrink=0.86)
        colorbar.set_label(r"minimum $\delta E$ (eV)")

    if float(log_velocity_ratio.min()) < 0 < float(log_velocity_ratio.max()):
        norm = TwoSlopeNorm(
            vmin=float(log_velocity_ratio.min()),
            vcenter=0.0,
            vmax=float(log_velocity_ratio.max()),
        )
    else:
        norm = None
    ratio_image = axes[2].imshow(
        log_velocity_ratio,
        origin="lower",
        extent=extent,
        aspect="equal",
        cmap="coolwarm",
        norm=norm,
        interpolation="nearest",
    )
    if float(log_velocity_ratio.min()) <= 0 <= float(log_velocity_ratio.max()):
        mu_mn_grid, mu_h2o_grid = np.meshgrid(mu_values, mu_values, indexing="xy")
        axes[2].contour(
            mu_mn_grid,
            mu_h2o_grid,
            log_velocity_ratio,
            levels=[0.0],
            colors="black",
            linewidths=1.8,
        )
    axes[2].set_title(
        rf"Growth kinetic ratio at {temperature:g} K" + "\n"
        + r"$\log_{10}(v_{001}/v_{011})$; black = equal"
    )
    axes[2].set_xlabel(r"$\Delta\mu_{\mathrm{Mn(OH)_2}}$ (eV)")
    axes[2].set_ylabel(r"$\Delta\mu_{\mathrm{H_2O}}$ (eV)")
    colorbar = fig.colorbar(ratio_image, ax=axes[2], shrink=0.86)
    colorbar.set_label(r"$\log_{10}(v_{001}/v_{011})$")
    if norm is not None:
        minimum = float(log_velocity_ratio.min())
        maximum = float(log_velocity_ratio.max())
        positive_ticks = np.linspace(0.0, maximum, 5)[1:]
        colorbar.set_ticks([minimum, 0.0, *positive_ticks])

    fig.suptitle(
        "Energetic-span-derived Mn$_3$O$_4$ growth comparison\n"
        r"$v=\Delta h(k_BT/h)\exp(-\delta E/k_BT)$; path boundaries from minimum-span sequence",
        fontsize=14,
    )
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.temperature <= 0:
        raise ValueError("temperature must be positive")
    if args.mu_step <= 0 or args.mu_max < args.mu_min:
        raise ValueError("invalid chemical-potential range")
    if args.height_001 <= 0 or args.height_011 <= 0:
        raise ValueError("repeat heights must be positive")

    spec_001 = SurfaceSpec(
        "001", args.energy_001, "S000", "S322", args.height_001
    )
    spec_011 = SurfaceSpec(
        "011", args.energy_011, "S000", "S644", args.height_011
    )
    prepared_001 = prepare_surface(spec_001)
    prepared_011 = prepare_surface(spec_011)

    count = int(round((args.mu_max - args.mu_min) / args.mu_step)) + 1
    mu_values = np.linspace(args.mu_min, args.mu_max, count)
    mu_mn_grid, mu_h2o_grid = np.meshgrid(mu_values, mu_values, indexing="xy")
    mu_mn_flat = mu_mn_grid.ravel()
    mu_h2o_flat = mu_h2o_grid.ravel()

    result_001 = evaluate_surface(
        prepared_001, mu_mn_flat, mu_h2o_flat, args.path_batch, args.grid_batch
    )
    result_011 = evaluate_surface(
        prepared_011, mu_mn_flat, mu_h2o_flat, args.path_batch, args.grid_batch
    )
    labels_001 = assign_region_labels(prepared_001, result_001)
    labels_011 = assign_region_labels(prepared_011, result_011)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    grid_csv = prefix.with_name(prefix.name + "_grid.csv")
    region_csv = prefix.with_name(prefix.name + "_path_regions.csv")
    boundary_csv = prefix.with_name(prefix.name + "_path_boundaries.csv")
    equal_growth_csv = prefix.with_name(prefix.name + "_equal_growth_boundaries.csv")
    figure_file = prefix.with_name(prefix.name + "_map.png")

    write_grid_csv(
        grid_csv,
        mu_mn_flat,
        mu_h2o_flat,
        prepared_001,
        result_001,
        labels_001,
        prepared_011,
        result_011,
        labels_011,
        args.temperature,
    )
    surfaces = [
        (prepared_001, result_001, labels_001),
        (prepared_011, result_011, labels_011),
    ]
    write_region_csv(region_csv, surfaces, mu_mn_flat, mu_h2o_flat)
    write_boundary_csv(
        boundary_csv, surfaces, mu_values, (len(mu_values), len(mu_values))
    )
    thermal_ev = KB_EV_PER_K * args.temperature
    log_ratio_grid = (
        math.log10(args.height_001 / args.height_011)
        + (result_011["best_span"] - result_001["best_span"])
        / (thermal_ev * math.log(10))
    ).reshape((len(mu_values), len(mu_values)))
    write_equal_growth_csv(equal_growth_csv, mu_values, log_ratio_grid)
    make_figure(
        figure_file,
        mu_values,
        prepared_001,
        result_001,
        labels_001,
        prepared_011,
        result_011,
        labels_011,
        args.temperature,
    )

    height_offset = math.log10(args.height_001 / args.height_011)
    log_ratio = (
        height_offset
        + (result_011["best_span"] - result_001["best_span"])
        / (thermal_ev * math.log(10))
    )
    print(f"Active [001] path regions: {len(labels_001)}")
    print(f"Active [011] path regions: {len(labels_011)}")
    print(
        "log10(v001/v011) range: "
        f"{log_ratio.min():.6g} to {log_ratio.max():.6g} at {args.temperature:g} K"
    )
    for output in (
        figure_file,
        grid_csv,
        region_csv,
        boundary_csv,
        equal_growth_csv,
    ):
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
