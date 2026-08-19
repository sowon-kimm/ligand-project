#!/usr/bin/env python3
"""Plot pure-to-capped states and subsequent capped-state reactions."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "book_capping.CSV"
DEFAULT_OUTPUT = SCRIPT_DIR / "pure_capping_connections.html"
DEFAULT_COMBINED_OUTPUT = SCRIPT_DIR / "pure_capping_connections_combined.html"
DEFAULT_COMBINED_MINIMUM_OUTPUT = (
    SCRIPT_DIR / "pure_capping_minimum_paths_combined.html"
)
CONTROL_TEMPLATE = SCRIPT_DIR / "pure_capping_controls.js"

REQUIRED_COLUMNS = {"state", "pure", "carboxyl", "amine"}
LIGANDS = ("carboxyl", "amine")
TRANSITIONS = ("OL", "C", "SR", "OX")
MU_LIMITS = (-4.0, 0.0)
LEVEL_HALF_WIDTH = 0.12
START_STATE = "S000"
FINAL_STATE = "S322"
TIE_TOLERANCE_EV = 1e-6

COLORS = {
    "pure": "#111111",
    "carboxyl": "#dc2626",
    "amine": "#2563eb",
    "OL": "#2563eb",
    "C": "#9333ea",
    "SR": "#16a34a",
    "OX": "#dc2626",
}


def parse_state(state: str) -> tuple[int, int, int]:
    """Parse S(O,C,X), requiring exactly three state counters."""
    if re.fullmatch(r"S\d{3}", state) is None:
        raise ValueError(f"Invalid state label {state!r}; expected a label such as S101")
    return tuple(int(digit) for digit in state[1:])


def reaction_coordinate(state: str) -> int:
    """Return the sum of the state counters, e.g. S101 -> 2."""
    return sum(parse_state(state))


def load_data(input_file: Path) -> pd.DataFrame:
    """Load and validate the pure/capping energy table."""
    if input_file.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(input_file)
    else:
        df = pd.read_csv(input_file)

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = df.copy()
    df["state"] = df["state"].astype(str).str.strip()
    counts = df["state"].map(parse_state)
    df[["O", "C", "X"]] = pd.DataFrame(counts.tolist(), index=df.index)
    df["reaction_coordinate"] = df[["O", "C", "X"]].sum(axis=1)
    df["pure_x"] = df["reaction_coordinate"].astype(float)
    df["capping_x"] = df["pure_x"] + 0.5

    for column in ("pure", "carboxyl", "amine"):
        df[column] = pd.to_numeric(df[column], errors="raise")

    if df["state"].duplicated().any():
        duplicates = sorted(df.loc[df["state"].duplicated(), "state"].unique())
        raise ValueError(f"Duplicate states found: {duplicates}")

    return df.sort_values(["reaction_coordinate", "state"]).reset_index(drop=True)


def classify_transition(a: dict, b: dict) -> str | None:
    """Classify an allowed S(O,C,X) counter change."""
    delta = (b["O"] - a["O"], b["C"] - a["C"], b["X"] - a["X"])
    return {
        (1, 0, 0): "OL",
        (0, 1, 0): "C",
        (0, 0, 1): "SR",
        (1, 1, 0): "OX",
    }.get(delta)


def build_edges(df: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Build the same reaction edges used by the ligand-project Plotly code."""
    states = df.to_dict("records")
    edges: list[tuple[str, str, str]] = []

    for a in states:
        for b in states:
            if b["reaction_coordinate"] <= a["reaction_coordinate"]:
                continue
            transition = classify_transition(a, b)
            if transition is None:
                continue
            expected_step = 2 if transition == "OX" else 1
            if b["reaction_coordinate"] - a["reaction_coordinate"] == expected_step:
                edges.append((a["state"], b["state"], transition))

    return edges


def enumerate_paths(
    edges: list[tuple[str, str, str]],
    start_state: str = START_STATE,
    final_state: str = FINAL_STATE,
) -> list[dict]:
    """Enumerate all complete reaction paths in stable lexical order."""
    adjacency: dict[str, list[tuple[str, str]]] = {}
    for source, target, transition in edges:
        adjacency.setdefault(source, []).append((target, transition))
    for source in adjacency:
        adjacency[source].sort()

    paths: list[dict] = []

    def dfs(state: str, states: list[str], mechanisms: list[str]) -> None:
        if state == final_state:
            paths.append(
                {"states": states.copy(), "mechanisms": mechanisms.copy()}
            )
            return
        for target, transition in adjacency.get(state, []):
            dfs(target, states + [target], mechanisms + [transition])

    dfs(start_state, [start_state], [])
    paths.sort(key=lambda path: " -> ".join(path["states"]))
    for index, path in enumerate(paths, start=1):
        path["path_id"] = f"P{index:03d}"
    if not paths:
        raise ValueError(f"No complete path from {start_state} to {final_state}")
    return paths


def corrected_energy(
    row: pd.Series | dict, energy_column: str, mu_mn: float, mu_h2o: float
) -> float:
    """Apply G = E0 - O*mu_Mn + (C + X/2)*mu_H2O."""
    return (
        row[energy_column]
        - row["O"] * mu_mn
        + (row["C"] + 0.5 * row["X"]) * mu_h2o
    )


def connection_arrays(
    df: pd.DataFrame, ligand: str, mu_mn: float = 0.0, mu_h2o: float = 0.0
) -> tuple[list, list]:
    """Build pure-to-capped line segments for one ligand."""
    x_values: list[float | None] = []
    y_values: list[float | None] = []

    for _, row in df.iterrows():
        x_values.extend(
            [
                row["pure_x"] + LEVEL_HALF_WIDTH,
                row["capping_x"] - LEVEL_HALF_WIDTH,
                None,
            ]
        )
        y_values.extend(
            [
                corrected_energy(row, "pure", mu_mn, mu_h2o),
                corrected_energy(row, ligand, mu_mn, mu_h2o),
                None,
            ]
        )

    return x_values, y_values


def reaction_edge_arrays(
    df: pd.DataFrame,
    edges: list[tuple[str, str, str]],
    ligand: str,
    transition: str,
    mu_mn: float = 0.0,
    mu_h2o: float = 0.0,
) -> tuple[list, list]:
    """Build capped-source to next-pure reaction segments."""
    state_map = df.set_index("state").to_dict("index")
    x_values: list[float | None] = []
    y_values: list[float | None] = []

    for from_name, to_name, edge_transition in edges:
        if edge_transition != transition:
            continue
        source = state_map[from_name]
        target = state_map[to_name]
        x_values.extend(
            [
                source["capping_x"] + LEVEL_HALF_WIDTH,
                target["pure_x"] - LEVEL_HALF_WIDTH,
                None,
            ]
        )
        y_values.extend(
            [
                corrected_energy(source, ligand, mu_mn, mu_h2o),
                corrected_energy(target, "pure", mu_mn, mu_h2o),
                None,
            ]
        )

    return x_values, y_values


def level_arrays(x_values, y_values) -> tuple[list, list]:
    """Build short horizontal energy levels."""
    level_x: list[float | None] = []
    level_y: list[float | None] = []
    for x_value, y_value in zip(x_values, y_values):
        level_x.extend(
            [x_value - LEVEL_HALF_WIDTH, x_value + LEVEL_HALF_WIDTH, None]
        )
        level_y.extend([y_value, y_value, None])
    return level_x, level_y


def corrected_series(
    df: pd.DataFrame, energy_column: str, mu_mn: float = 0.0, mu_h2o: float = 0.0
) -> list[float]:
    return [
        corrected_energy(row, energy_column, mu_mn, mu_h2o)
        for _, row in df.iterrows()
    ]


def add_trace(fig: go.Figure, trace: go.Scatter, row: int) -> int:
    """Add a trace to a subplot and return its global trace index."""
    trace_index = len(fig.data)
    fig.add_trace(trace, row=row, col=1)
    return trace_index


def make_figure(
    df: pd.DataFrame, edges: list[tuple[str, str, str]]
) -> tuple[go.Figure, dict]:
    """Create two ligand panels and return trace indices for JS updates."""
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0.09,
        subplot_titles=("Carboxyl capping", "Amine capping"),
    )
    trace_indices: dict[str, dict] = {}

    for row_number, ligand in enumerate(LIGANDS, start=1):
        trace_indices[ligand] = {"reaction": {}}
        link_x, link_y = connection_arrays(df, ligand)
        trace_indices[ligand]["capping_link"] = add_trace(
            fig,
            go.Scatter(
                x=link_x,
                y=link_y,
                mode="lines",
                line={
                    "width": 2.0,
                    "color": COLORS["pure"],
                    "dash": "solid" if ligand == "carboxyl" else "dash",
                },
                opacity=0.65,
                name="Pure → capped",
                legendgroup="capping_link",
                showlegend=row_number == 1,
                hoverinfo="skip",
            ),
            row_number,
        )

        for transition in TRANSITIONS:
            edge_x, edge_y = reaction_edge_arrays(
                df, edges, ligand, transition
            )
            trace_indices[ligand]["reaction"][transition] = add_trace(
                fig,
                go.Scatter(
                    x=edge_x,
                    y=edge_y,
                    mode="lines",
                    line={"width": 2.2, "color": COLORS[transition]},
                    opacity=0.72,
                    name=transition,
                    legendgroup=f"reaction_{transition}",
                    showlegend=row_number == 1,
                    hoverinfo="skip",
                ),
                row_number,
            )

        pure_y = corrected_series(df, "pure")
        capped_y = corrected_series(df, ligand)
        pure_level_x, pure_level_y = level_arrays(df["pure_x"], pure_y)
        capped_level_x, capped_level_y = level_arrays(df["capping_x"], capped_y)

        trace_indices[ligand]["pure_levels"] = add_trace(
            fig,
            go.Scatter(
                x=pure_level_x,
                y=pure_level_y,
                mode="lines",
                line={"width": 3.5, "color": COLORS["pure"]},
                hoverinfo="skip",
                showlegend=False,
            ),
            row_number,
        )
        trace_indices[ligand]["capped_levels"] = add_trace(
            fig,
            go.Scatter(
                x=capped_level_x,
                y=capped_level_y,
                mode="lines",
                line={"width": 3.5, "color": COLORS[ligand]},
                hoverinfo="skip",
                showlegend=False,
            ),
            row_number,
        )

        pure_custom = df[["state", "reaction_coordinate", "pure", "O", "C", "X"]]
        trace_indices[ligand]["pure_markers"] = add_trace(
            fig,
            go.Scatter(
                x=df["pure_x"],
                y=pure_y,
                mode="markers+text",
                text=df["state"],
                textposition="top center",
                textfont={"size": 12, "color": COLORS["pure"]},
                marker={
                    "size": 7,
                    "symbol": "circle",
                    "color": COLORS["pure"],
                    "line": {"width": 0.8, "color": "white"},
                },
                name="Pure state",
                legendgroup="pure",
                showlegend=row_number == 1,
                customdata=pure_custom.to_numpy(),
                hovertemplate=(
                    "<b>%{customdata[0]} — pure</b><br>"
                    "Reaction coordinate: %{customdata[1]}<br>"
                    "Corrected energy: %{y:.6f} eV<br>"
                    "E0: %{customdata[2]:.6f} eV<br>"
                    "O=%{customdata[3]}, C=%{customdata[4]}, X=%{customdata[5]}"
                    "<extra></extra>"
                ),
            ),
            row_number,
        )

        capped_custom = pd.DataFrame(
            {
                "state": df["state"],
                "coordinate": df["reaction_coordinate"],
                "E0": df[ligand],
                "delta": df[ligand] - df["pure"],
                "O": df["O"],
                "C": df["C"],
                "X": df["X"],
            }
        )
        trace_indices[ligand]["capped_markers"] = add_trace(
            fig,
            go.Scatter(
                x=df["capping_x"],
                y=capped_y,
                mode="markers",
                marker={
                    "size": 9,
                    "symbol": "diamond" if ligand == "carboxyl" else "triangle-up",
                    "color": COLORS[ligand],
                    "line": {"width": 0.8, "color": "white"},
                },
                name="Ligand-capped state",
                legendgroup="capped",
                showlegend=row_number == 1,
                customdata=capped_custom.to_numpy(),
                hovertemplate=(
                    f"<b>%{{customdata[0]}} — {ligand}</b><br>"
                    "Reaction coordinate: %{customdata[1]}<br>"
                    "Corrected energy: %{y:.6f} eV<br>"
                    "E0: %{customdata[2]:.6f} eV<br>"
                    "ΔE(capping − pure): %{customdata[3]:.6f} eV<br>"
                    "O=%{customdata[4]}, C=%{customdata[5]}, X=%{customdata[6]}"
                    "<extra></extra>"
                ),
            ),
            row_number,
        )

    fixed_energies = [
        corrected_energy(row, energy_column, mu_mn, mu_h2o)
        for _, row in df.iterrows()
        for energy_column in ("pure", "carboxyl", "amine")
        for mu_mn in MU_LIMITS
        for mu_h2o in MU_LIMITS
    ]
    energy_padding = max(0.5, 0.05 * (max(fixed_energies) - min(fixed_energies)))
    y_range = [min(fixed_energies) - energy_padding, max(fixed_energies) + energy_padding]
    max_coordinate = int(df["reaction_coordinate"].max())

    fig.update_layout(
        title=(
            "<b>Chemical-potential dependent ligand-capped reaction networks</b>"
            "<br><sup>ΔμMn(OH)₂ = 0.00 eV, ΔμH₂O = 0.00 eV</sup>"
        ),
        template="plotly_white",
        autosize=False,
        width=1400,
        height=1250,
        margin={"l": 90, "r": 45, "t": 150, "b": 80},
        hovermode="closest",
        legend={
            "orientation": "h",
            "x": 0.58,
            "xanchor": "center",
            "y": 1.08,
        },
    )
    fig.update_xaxes(
        title_text="Reaction coordinate",
        dtick=0.5,
        range=[-0.35, max_coordinate + 0.85],
        fixedrange=False,
        row=2,
        col=1,
    )
    fig.update_xaxes(
        dtick=0.5,
        range=[-0.35, max_coordinate + 0.85],
        fixedrange=False,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title_text="Relative energy (eV)",
        range=y_range,
        autorange=False,
        fixedrange=False,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title_text="Relative energy (eV)",
        range=y_range,
        autorange=False,
        fixedrange=False,
        row=2,
        col=1,
    )

    return fig, trace_indices


def make_combined_figure(
    df: pd.DataFrame, edges: list[tuple[str, str, str]]
) -> tuple[go.Figure, dict]:
    """Overlay the carboxyl and amine capped networks in one coordinate system."""
    fig = make_subplots(rows=1, cols=1)
    trace_indices: dict[str, dict] = {}

    for ligand in LIGANDS:
        ligand_label = ligand.capitalize()
        trace_indices[ligand] = {"reaction": {}}
        link_x, link_y = connection_arrays(df, ligand)
        trace_indices[ligand]["capping_link"] = add_trace(
            fig,
            go.Scatter(
                x=link_x,
                y=link_y,
                mode="lines",
                line={
                    "width": 2.0,
                    "color": COLORS["pure"],
                    "dash": "solid" if ligand == "carboxyl" else "dash",
                },
                opacity=0.65,
                name=f"{ligand_label} adsorption",
                legendgroup=f"combined_{ligand}",
                hoverinfo="skip",
            ),
            1,
        )

        for transition in TRANSITIONS:
            edge_x, edge_y = reaction_edge_arrays(
                df, edges, ligand, transition
            )
            trace_indices[ligand]["reaction"][transition] = add_trace(
                fig,
                go.Scatter(
                    x=edge_x,
                    y=edge_y,
                    mode="lines",
                    line={
                        "width": 2.2,
                        "color": COLORS[transition],
                        "dash": "solid" if ligand == "carboxyl" else "dash",
                    },
                    opacity=0.72,
                    name=f"{ligand_label} · {transition}",
                    legendgroup=f"combined_{ligand}",
                    hoverinfo="skip",
                ),
                1,
            )

        capped_y = corrected_series(df, ligand)
        capped_level_x, capped_level_y = level_arrays(df["capping_x"], capped_y)
        trace_indices[ligand]["capped_levels"] = add_trace(
            fig,
            go.Scatter(
                x=capped_level_x,
                y=capped_level_y,
                mode="lines",
                line={"width": 3.5, "color": COLORS[ligand]},
                hoverinfo="skip",
                showlegend=False,
            ),
            1,
        )

        capped_custom = pd.DataFrame(
            {
                "state": df["state"],
                "coordinate": df["reaction_coordinate"],
                "E0": df[ligand],
                "delta": df[ligand] - df["pure"],
                "O": df["O"],
                "C": df["C"],
                "X": df["X"],
            }
        )
        trace_indices[ligand]["capped_markers"] = add_trace(
            fig,
            go.Scatter(
                x=df["capping_x"],
                y=capped_y,
                mode="markers",
                marker={
                    "size": 9,
                    "symbol": "diamond" if ligand == "carboxyl" else "triangle-up",
                    "color": COLORS[ligand],
                    "line": {"width": 0.8, "color": "white"},
                },
                name=f"{ligand_label} capped state",
                legendgroup=f"combined_{ligand}",
                showlegend=False,
                customdata=capped_custom.to_numpy(),
                hovertemplate=(
                    f"<b>%{{customdata[0]}} — {ligand}</b><br>"
                    "Reaction coordinate: %{customdata[1]}<br>"
                    "Corrected energy: %{y:.6f} eV<br>"
                    "E0: %{customdata[2]:.6f} eV<br>"
                    "ΔE(capping − pure): %{customdata[3]:.6f} eV<br>"
                    "O=%{customdata[4]}, C=%{customdata[5]}, X=%{customdata[6]}"
                    "<extra></extra>"
                ),
            ),
            1,
        )

    pure_y = corrected_series(df, "pure")
    pure_level_x, pure_level_y = level_arrays(df["pure_x"], pure_y)
    pure_level_index = add_trace(
        fig,
        go.Scatter(
            x=pure_level_x,
            y=pure_level_y,
            mode="lines",
            line={"width": 3.5, "color": COLORS["pure"]},
            hoverinfo="skip",
            showlegend=False,
        ),
        1,
    )
    pure_custom = df[["state", "reaction_coordinate", "pure", "O", "C", "X"]]
    pure_marker_index = add_trace(
        fig,
        go.Scatter(
            x=df["pure_x"],
            y=pure_y,
            mode="markers+text",
            text=df["state"],
            textposition="top center",
            textfont={"size": 12, "color": COLORS["pure"]},
            marker={
                "size": 7,
                "symbol": "circle",
                "color": COLORS["pure"],
                "line": {"width": 0.8, "color": "white"},
            },
            name="Pure state",
            legendgroup="pure",
            customdata=pure_custom.to_numpy(),
            hovertemplate=(
                "<b>%{customdata[0]} — pure</b><br>"
                "Reaction coordinate: %{customdata[1]}<br>"
                "Corrected energy: %{y:.6f} eV<br>"
                "E0: %{customdata[2]:.6f} eV<br>"
                "O=%{customdata[3]}, C=%{customdata[4]}, X=%{customdata[5]}"
                "<extra></extra>"
            ),
        ),
        1,
    )
    for ligand in LIGANDS:
        trace_indices[ligand]["pure_levels"] = pure_level_index
        trace_indices[ligand]["pure_markers"] = pure_marker_index

    fixed_energies = [
        corrected_energy(row, energy_column, mu_mn, mu_h2o)
        for _, row in df.iterrows()
        for energy_column in ("pure", "carboxyl", "amine")
        for mu_mn in MU_LIMITS
        for mu_h2o in MU_LIMITS
    ]
    energy_padding = max(0.5, 0.05 * (max(fixed_energies) - min(fixed_energies)))
    max_coordinate = int(df["reaction_coordinate"].max())

    fig.update_layout(
        title=(
            "<b>Combined carboxyl and amine capped reaction networks</b>"
            "<br><sup>Solid reaction lines: carboxyl | dashed reaction lines: amine"
            " | ΔμMn(OH)₂ = 0.00 eV, ΔμH₂O = 0.00 eV</sup>"
        ),
        template="plotly_white",
        autosize=False,
        width=1400,
        height=900,
        margin={"l": 90, "r": 45, "t": 155, "b": 80},
        hovermode="closest",
        legend={
            "orientation": "h",
            "x": 0.58,
            "xanchor": "center",
            "y": 1.15,
        },
    )
    fig.update_xaxes(
        title_text="Reaction coordinate",
        dtick=0.5,
        range=[-0.35, max_coordinate + 0.85],
        fixedrange=False,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title_text="Relative energy (eV)",
        range=[
            min(fixed_energies) - energy_padding,
            max(fixed_energies) + energy_padding,
        ],
        autorange=False,
        fixedrange=False,
        row=1,
        col=1,
    )

    return fig, trace_indices


def make_combined_minimum_figure(
    df: pd.DataFrame, edges: list[tuple[str, str, str]]
) -> tuple[go.Figure, dict]:
    """Create the combined canvas used for dynamically selected minimum paths."""
    fig, trace_indices = make_combined_figure(df, edges)
    trace_indices["carboxyl"]["tdi_marker"] = add_trace(
        fig,
        go.Scatter(
            x=[], y=[], mode="markers+text",
            marker={
                "size": 17,
                "symbol": "diamond",
                "color": "#00AEEF",
                "line": {"width": 2, "color": COLORS["carboxyl"]},
            },
            textposition="bottom center",
            textfont={"size": 11, "color": COLORS["carboxyl"]},
            name="Carboxyl TDI",
            zorder=100,
            cliponaxis=False,
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>TDI for carboxyl<br>"
                "Corrected energy: %{y:.6f} eV<extra></extra>"
            ),
        ),
        1,
    )
    trace_indices["carboxyl"]["tdts_marker"] = add_trace(
        fig,
        go.Scatter(
            x=[], y=[], mode="markers+text",
            marker={
                "size": 21,
                "symbol": "star",
                "color": "#E53935",
                "line": {"width": 2, "color": COLORS["carboxyl"]},
            },
            textposition="bottom center",
            textfont={"size": 11, "color": COLORS["carboxyl"]},
            name="Carboxyl TDTS proxy",
            zorder=101,
            cliponaxis=False,
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>TDTS proxy for carboxyl<br>"
                "Corrected energy: %{y:.6f} eV<extra></extra>"
            ),
        ),
        1,
    )
    fig.update_layout(
        title=(
            "<b>Minimum pure → capped half-step → next pure pathways</b>"
            "<br><sup>TDTS is an energy-level proxy</sup>"
        ),
        height=950,
        margin={"l": 90, "r": 45, "t": 205, "b": 80},
        legend={
            "orientation": "h",
            "x": 0.58,
            "xanchor": "center",
            "y": 1.20,
        },
    )
    return fig, trace_indices


def state_payload(df: pd.DataFrame) -> list[dict]:
    columns = [
        "state",
        "reaction_coordinate",
        "pure_x",
        "capping_x",
        "pure",
        "carboxyl",
        "amine",
        "O",
        "C",
        "X",
    ]
    states = df[columns].rename(columns={"reaction_coordinate": "coordinate"})
    return states.to_dict("records")


def write_interactive_html(
    fig: go.Figure,
    output_file: Path,
    df: pd.DataFrame,
    edges: list[tuple[str, str, str]],
    trace_indices: dict,
    layout_mode: str,
    all_paths: list[dict],
) -> None:
    """Write a standalone Plotly HTML file with two chemical-potential sliders."""
    payload = {
        "states": state_payload(df),
        "edges": edges,
        "ligands": list(LIGANDS),
        "transitions": list(TRANSITIONS),
        "levelHalfWidth": LEVEL_HALF_WIDTH,
        "traceIndices": trace_indices,
        "layoutMode": layout_mode,
        "allPaths": all_paths,
        "tieTolerance": TIE_TOLERANCE_EV,
        "muMin": MU_LIMITS[0],
        "muMax": MU_LIMITS[1],
    }
    post_script = CONTROL_TEMPLATE.read_text(encoding="utf-8").replace(
        "__PAYLOAD_JSON__",
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
    )
    fig.write_html(
        output_file,
        include_plotlyjs=True,
        full_html=True,
        config={
            "responsive": False,
            "scrollZoom": True,
            "displayModeBar": True,
            "displaylogo": False,
            "doubleClick": "reset",
        },
        post_script=post_script,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input CSV or Excel file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--layout",
        choices=("panels", "combined", "combined-minimum"),
        default="panels",
        help=(
            "Use separate panels, overlay all paths, or overlay only minimum "
            "energetic-span paths (default: panels)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output standalone HTML file (default depends on --layout)",
    )
    parser.add_argument(
        "--start-state",
        default=START_STATE,
        help=f"Initial state for complete paths (default: {START_STATE})",
    )
    parser.add_argument(
        "--final-state",
        default=FINAL_STATE,
        help=f"Final state for complete paths (default: {FINAL_STATE})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = load_data(args.input)
    edges = build_edges(df)
    all_paths = enumerate_paths(edges, args.start_state, args.final_state)
    if args.layout == "combined-minimum":
        fig, trace_indices = make_combined_minimum_figure(df, edges)
        output_file = args.output or DEFAULT_COMBINED_MINIMUM_OUTPUT
    elif args.layout == "combined":
        fig, trace_indices = make_combined_figure(df, edges)
        output_file = args.output or DEFAULT_COMBINED_OUTPUT
    else:
        fig, trace_indices = make_figure(df, edges)
        output_file = args.output or DEFAULT_OUTPUT
    output_file.parent.mkdir(parents=True, exist_ok=True)
    write_interactive_html(
        fig, output_file, df, edges, trace_indices, args.layout, all_paths
    )
    print(f"Loaded {len(df)} states from: {args.input}")
    print(f"Built {len(edges)} subsequent reaction edges per ligand")
    print(f"Enumerated {len(all_paths)} complete paths")
    print(f"Layout: {args.layout}")
    print(f"Saved interactive figure to: {output_file}")


if __name__ == "__main__":
    main()
