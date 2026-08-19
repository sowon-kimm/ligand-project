#!/usr/bin/env python3
"""Plot minimum paths with pure integer and capped half-step coordinates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import plotly.graph_objects as go

from plot_pure_capping_connections import (
    COLORS,
    LIGANDS,
    MU_LIMITS,
    TRANSITIONS,
    build_edges,
    enumerate_paths,
    load_data,
    state_payload,
)


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "book_capping.CSV"
DEFAULT_OUTPUT = HERE / "pure_capping_desorb_react_adsorb_minimum_001.html"
CONTROL_TEMPLATE = HERE / "desorb_react_adsorb_controls.js"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start-state", default="S000")
    parser.add_argument("--final-state", default="S322")
    return parser.parse_args()


def add_trace(fig: go.Figure, trace: go.Scatter) -> int:
    index = len(fig.data)
    fig.add_trace(trace)
    return index


def make_figure(df) -> tuple[go.Figure, dict]:
    fig = go.Figure()
    indices: dict[str, dict] = {}

    for ligand in LIGANDS:
        label = ligand.capitalize()
        ligand_dash = "solid" if ligand == "carboxyl" else "dash"
        indices[ligand] = {"reaction": {}}
        indices[ligand]["desorption"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="lines",
                line={
                    "width": 2.5,
                    "color": "#111111",
                    "dash": "solid" if ligand == "carboxyl" else "dash",
                },
                name=f"{label} adsorption (half-step)", legendgroup=ligand,
                hoverinfo="skip",
            ),
        )
        for transition in TRANSITIONS:
            indices[ligand]["reaction"][transition] = add_trace(
                fig,
                go.Scatter(
                    x=[], y=[], mode="lines",
                    line={
                        "width": 3.2,
                        "color": COLORS[transition],
                        "dash": ligand_dash,
                    },
                    name=f"{label} · {transition} (desorb + react)",
                    legendgroup=ligand,
                    hoverinfo="skip",
                ),
            )
        indices[ligand]["adsorption"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="lines",
                line={"width": 2.5, "color": COLORS[ligand], "dash": "dot"},
                showlegend=False, visible=False, hoverinfo="skip",
            ),
        )
        indices[ligand]["capped_levels"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="lines",
                line={"width": 4, "color": COLORS[ligand]},
                showlegend=False, hoverinfo="skip",
            ),
        )
        indices[ligand]["capped_markers"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="markers",
                marker={
                    "size": 10,
                    "symbol": "diamond" if ligand == "carboxyl" else "triangle-up",
                    "color": COLORS[ligand],
                    "line": {"width": 0.8, "color": "white"},
                },
                showlegend=False,
                hovertemplate=(
                    "<b>%{customdata[0]} — %{customdata[1]}</b><br>"
                    "Path step: %{x:.2f}<br>"
                    "Corrected energy: %{y:.6f} eV<br>"
                    "Step ΔE: %{customdata[2]:.6f} eV"
                    "<extra></extra>"
                ),
            ),
        )
        indices[ligand]["pure_levels"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="lines",
                line={"width": 3, "color": "#111111", "dash": ligand_dash},
                opacity=0.75, showlegend=False, hoverinfo="skip",
            ),
        )
        indices[ligand]["pure_markers"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="markers+text",
                marker={
                    "size": 8,
                    "symbol": "circle",
                    "color": "#111111",
                    "line": {"width": 0.8, "color": "white"},
                },
                textposition="top center",
                textfont={"size": 12, "color": "#111111"},
                showlegend=False,
                hovertemplate=(
                    "<b>%{customdata[0]} — %{customdata[1]}</b><br>"
                    "Path step: %{x:.2f}<br>"
                    "Corrected pure energy: %{y:.6f} eV<br>"
                    "Step ΔE: %{customdata[2]:.6f} eV"
                    "<extra></extra>"
                ),
            ),
        )
        indices[ligand]["tdi_marker"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="markers+text",
                marker={
                    "size": 17,
                    "symbol": "diamond",
                    "color": "#00AEEF",
                    "line": {"width": 2, "color": COLORS[ligand]},
                },
                textposition="bottom center",
                textfont={"size": 11, "color": COLORS[ligand]},
                name=f"{label} TDI",
                legendgroup=f"{ligand}_span",
                visible=ligand == "carboxyl",
                showlegend=ligand == "carboxyl",
                zorder=100,
                cliponaxis=False,
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "TDI for %{customdata[1]}<br>"
                    "Corrected energy: %{y:.6f} eV"
                    "<extra></extra>"
                ),
            ),
        )
        indices[ligand]["tdts_marker"] = add_trace(
            fig,
            go.Scatter(
                x=[], y=[], mode="markers+text",
                marker={
                    "size": 21,
                    "symbol": "star",
                    "color": "#E53935",
                    "line": {"width": 2, "color": COLORS[ligand]},
                },
                textposition="bottom center",
                textfont={"size": 11, "color": COLORS[ligand]},
                name=f"{label} TDTS proxy",
                legendgroup=f"{ligand}_span",
                visible=ligand == "carboxyl",
                showlegend=ligand == "carboxyl",
                zorder=101,
                cliponaxis=False,
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "TDTS proxy for %{customdata[1]}<br>"
                    "Corrected energy: %{y:.6f} eV"
                    "<extra></extra>"
                ),
            ),
        )

    fixed = []
    for row in df.to_dict("records"):
        for key in ("pure", "carboxyl", "amine"):
            for mu_mn in MU_LIMITS:
                for mu_h2o in MU_LIMITS:
                    fixed.append(
                        row[key]
                        - row["O"] * mu_mn
                        + (row["C"] + 0.5 * row["X"]) * mu_h2o
                    )
    padding = max(0.5, 0.05 * (max(fixed) - min(fixed)))
    fig.update_layout(
        title=(
            "<b>Pure integer → capped half-step → next pure pathways (001)</b>"
            "<br><sup>Desorption–reaction–adsorption model; TDTS is an energy-level proxy</sup>"
        ),
        template="plotly_white",
        autosize=False,
        width=1500,
        height=920,
        margin={"l": 90, "r": 45, "t": 220, "b": 90},
        hovermode="closest",
        legend={
            "orientation": "h", "x": 0.58, "xanchor": "center", "y": 1.22,
        },
        xaxis={
            "title": "Reaction coordinate (pure integer → capped half-step → next pure integer)",
            "dtick": 0.5,
            "range": [-0.35, 8.35],
            "fixedrange": False,
        },
        yaxis={
            "title": "Relative energy (eV)",
            "range": [min(fixed) - padding, max(fixed) + padding],
            "autorange": False,
            "fixedrange": False,
        },
    )
    return fig, indices


def main() -> None:
    args = parse_args()
    df = load_data(args.input)
    edges = build_edges(df)
    paths = enumerate_paths(edges, args.start_state, args.final_state)
    fig, trace_indices = make_figure(df)
    payload = {
        "states": state_payload(df),
        "paths": paths,
        "ligands": list(LIGANDS),
        "transitions": list(TRANSITIONS),
        "traceIndices": trace_indices,
        "levelHalfWidth": 0.09,
        "tieTolerance": 1e-6,
        "muMin": MU_LIMITS[0],
        "muMax": MU_LIMITS[1],
    }
    post_script = CONTROL_TEMPLATE.read_text(encoding="utf-8").replace(
        "__PAYLOAD_JSON__",
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        args.output,
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
    print(f"Loaded {len(df)} states and enumerated {len(paths)} complete paths")
    print(f"Saved desorption/reaction/adsorption minimum-path figure to: {args.output}")


if __name__ == "__main__":
    main()
