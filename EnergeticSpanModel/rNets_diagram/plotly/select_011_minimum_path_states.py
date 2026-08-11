#!/usr/bin/env python3
"""Select all 011 layer-1--4 states on minimum energetic-span paths."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_ENERGY = HERE / "energy_011_1234.csv"
DEFAULT_OUTPUT = HERE / "minimum_path_states_011_1234_mu_-4_0.csv"
START_STATE = "S000"
FINAL_STATE = "S644"
TIE_TOLERANCE_EV = 1e-6


def load_states(path: Path) -> list[dict]:
    states = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            state = row["state"].strip()
            energy = row["relative_E"].strip()
            if not state or not energy:
                continue
            if not re_state(state):
                raise ValueError(f"Invalid state label: {state}")
            o_count, c_count, x_count = map(int, state[1:])
            states.append(
                {
                    "path": row.get("path", "").strip(),
                    "state": state,
                    "coordinate": float(row["coordinate"]),
                    "E0": float(energy),
                    "O": o_count,
                    "C": c_count,
                    "X": x_count,
                }
            )
    return states


def re_state(value: str) -> bool:
    return len(value) == 4 and value.startswith("S") and value[1:].isdigit()


def transition(a: dict, b: dict) -> str | None:
    delta = (b["O"] - a["O"], b["C"] - a["C"], b["X"] - a["X"])
    return {(1, 0, 0): "OL", (0, 1, 0): "C", (0, 0, 1): "SR", (1, 1, 0): "OX"}.get(delta)


def build_edges(states: list[dict]) -> list[tuple[str, str, str]]:
    edges = []
    for a in states:
        for b in states:
            if b["coordinate"] <= a["coordinate"]:
                continue
            mechanism = transition(a, b)
            if mechanism is None:
                continue
            expected = 2 if mechanism == "OX" else 1
            if abs((b["coordinate"] - a["coordinate"]) - expected) < 1e-9:
                edges.append((a["state"], b["state"], mechanism))
    return edges


def enumerate_paths(
    edges: list[tuple[str, str, str]], start_state: str, final_state: str
) -> list[list[str]]:
    adjacency: dict[str, list[str]] = defaultdict(list)
    for start, end, _ in edges:
        adjacency[start].append(end)
    for node in adjacency:
        adjacency[node].sort()
    paths = []

    def dfs(node: str, current: list[str]) -> None:
        if node == final_state:
            paths.append(current.copy())
            return
        for next_node in adjacency.get(node, []):
            if next_node not in current:
                dfs(next_node, current + [next_node])

    dfs(start_state, [start_state])
    paths.sort(key=lambda path: " -> ".join(path))
    return paths


def path_spans(
    candidate_indices: np.ndarray,
    candidate_valid: np.ndarray,
    state_energies: np.ndarray,
    grid_indices: np.ndarray,
    reaction_energy: np.ndarray,
) -> np.ndarray:
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
        rise = energy - new_minimum
        drop = maximum_prefix - energy
        maximum_rise = np.where(valid, np.maximum(maximum_rise, rise), maximum_rise)
        maximum_drop = np.where(valid, np.maximum(maximum_drop, drop), maximum_drop)
        minimum_prefix = np.where(valid, new_minimum, minimum_prefix)
        maximum_prefix = np.where(valid, np.maximum(maximum_prefix, energy), maximum_prefix)
    return np.maximum(maximum_rise, maximum_drop + reaction_energy[grid_indices][None, :])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy", type=Path, default=DEFAULT_ENERGY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mu-min", type=float, default=-4.0)
    parser.add_argument("--mu-max", type=float, default=0.0)
    parser.add_argument("--mu-step", type=float, default=0.05)
    parser.add_argument("--start-state", default=START_STATE)
    parser.add_argument("--final-state", default=FINAL_STATE)
    parser.add_argument("--path-batch", type=int, default=4000)
    parser.add_argument("--grid-batch", type=int, default=243)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Output exists; refusing to overwrite: {args.output}")

    states = load_states(args.energy)
    state_index = {state["state"]: index for index, state in enumerate(states)}
    if args.start_state not in state_index or args.final_state not in state_index:
        raise ValueError("Start or final state missing")
    edges = build_edges(states)
    paths = enumerate_paths(edges, args.start_state, args.final_state)
    if not paths:
        raise ValueError("No complete paths")
    maximum_candidates = max(len(path) - 1 for path in paths)
    candidate_indices = np.full((len(paths), maximum_candidates), -1, dtype=np.int16)
    path_state_indices = []
    for row, path in enumerate(paths):
        indices = [state_index[name] for name in path]
        path_state_indices.append(np.array(indices, dtype=np.int16))
        candidate_indices[row, : len(indices) - 1] = indices[:-1]
    candidate_valid = candidate_indices >= 0

    count = int(round((args.mu_max - args.mu_min) / args.mu_step)) + 1
    mu_values = np.linspace(args.mu_min, args.mu_max, count)
    mu_mn, mu_h2o = np.meshgrid(mu_values, mu_values, indexing="ij")
    mu_mn = mu_mn.ravel()
    mu_h2o = mu_h2o.ravel()
    state_energies = np.array(
        [
            state["E0"] - state["O"] * mu_mn + (state["C"] + 0.5 * state["X"]) * mu_h2o
            for state in states
        ]
    )
    reaction_energy = (
        state_energies[state_index[args.final_state]]
        - state_energies[state_index[args.start_state]]
    )
    grid_count = len(mu_mn)
    best_span = np.full(grid_count, np.inf)
    representative_path = np.full(grid_count, -1, dtype=np.int32)

    for grid_start in range(0, grid_count, args.grid_batch):
        grid_indices = np.arange(grid_start, min(grid_start + args.grid_batch, grid_count))
        local_best = np.full(len(grid_indices), np.inf)
        local_path = np.full(len(grid_indices), -1, dtype=np.int32)
        for path_start in range(0, len(paths), args.path_batch):
            path_end = min(path_start + args.path_batch, len(paths))
            spans = path_spans(
                candidate_indices[path_start:path_end], candidate_valid[path_start:path_end],
                state_energies, grid_indices, reaction_energy,
            )
            indices = np.argmin(spans, axis=0)
            values = spans[indices, np.arange(len(grid_indices))]
            better = values < local_best - TIE_TOLERANCE_EV
            local_best[better] = values[better]
            local_path[better] = path_start + indices[better]
        best_span[grid_indices] = local_best
        representative_path[grid_indices] = local_path
        print(f"minimum pass: {grid_indices[-1] + 1}/{grid_count}", flush=True)

    cooptimal_on_grid = np.zeros((len(states), grid_count), dtype=bool)
    for grid_start in range(0, grid_count, args.grid_batch):
        grid_indices = np.arange(grid_start, min(grid_start + args.grid_batch, grid_count))
        for path_start in range(0, len(paths), args.path_batch):
            path_end = min(path_start + args.path_batch, len(paths))
            spans = path_spans(
                candidate_indices[path_start:path_end], candidate_valid[path_start:path_end],
                state_energies, grid_indices, reaction_energy,
            )
            tied = np.abs(spans - best_span[grid_indices][None, :]) <= TIE_TOLERANCE_EV
            for local_path in np.flatnonzero(tied.any(axis=1)):
                selected_grids = grid_indices[tied[local_path]]
                selected_states = path_state_indices[path_start + local_path]
                cooptimal_on_grid[np.ix_(selected_states, selected_grids)] = True
        print(f"cooptimal pass: {grid_indices[-1] + 1}/{grid_count}", flush=True)

    representative_on_grid = np.zeros((len(states), grid_count), dtype=bool)
    for grid, path_index in enumerate(representative_path):
        representative_on_grid[path_state_indices[path_index], grid] = True

    rows = []
    for index, state in enumerate(states):
        cooptimal_count = int(cooptimal_on_grid[index].sum())
        representative_count = int(representative_on_grid[index].sum())
        if cooptimal_count == 0:
            continue
        rows.append(
            {
                "state": state["state"],
                "path": state["path"],
                "coordinate": f"{state['coordinate']:g}",
                "relative_E": f"{state['E0']:.9f}",
                "cooptimal_grid_points": cooptimal_count,
                "representative_grid_points": representative_count,
                "selection": "representative" if representative_count else "cooptimal_only",
            }
        )
    columns = (
        "state", "path", "coordinate", "relative_E", "cooptimal_grid_points",
        "representative_grid_points", "selection",
    )
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Selected {len(rows)} states from {len(paths)} complete paths on "
        f"{grid_count} chemical-potential points -> {args.output}"
    )


if __name__ == "__main__":
    main()
