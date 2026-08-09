from pathlib import Path
import csv

import plotly.graph_objects as go
import ipywidgets as widgets
from IPython.display import display


# ============================================================
# 1. Input CSV
# ============================================================
# Put energy_001.csv in the same directory as this script.
CSV_FILE = Path(__file__).resolve().with_name("energy_001.csv")


def load_states(csv_file):
    """
    Read:
        path,state,coordinate,relative_E

    O, C, X are parsed automatically from state labels such as S211:
        S211 -> O=2, C=1, X=1

    Rows with blank relative_E are skipped.
    """
    data = []

    with open(csv_file, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        required = {"state", "coordinate", "relative_E"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing required CSV columns: {sorted(missing)}\n"
                f"Required columns: state, coordinate, relative_E"
            )

        for row in reader:
            state = row["state"].strip()
            energy_text = row["relative_E"].strip()

            # Allow incomplete states, e.g. S002 with no energy yet.
            if not state or not energy_text:
                continue

            if len(state) != 4 or not state.startswith("S") or not state[1:].isdigit():
                raise ValueError(
                    f"State '{state}' must have the form S(O,C,X), e.g. S211."
                )

            O, C, X = map(int, state[1:])

            data.append(
                {
                    "path": row.get("path", "").strip(),
                    "state": state,
                    "coordinate": float(row["coordinate"]),
                    "E0": float(energy_text),
                    "O": O,
                    "C": C,
                    "X": X,
                }
            )

    return data


data = load_states(CSV_FILE)
smap = {s["state"]: s for s in data}


# ============================================================
# 2. Reaction network
# ============================================================
# OL : O + 1
# C  : C + 1
# SR : X + 1
# OX : O + 1 and C + 1
#
# Edges are generated automatically from the loaded states,
# so you do not need to manually list every transition.
# ============================================================
def classify_transition(a, b):
    dO = b["O"] - a["O"]
    dC = b["C"] - a["C"]
    dX = b["X"] - a["X"]

    if (dO, dC, dX) == (1, 0, 0):
        return "OL"
    if (dO, dC, dX) == (0, 1, 0):
        return "C"
    if (dO, dC, dX) == (0, 0, 1):
        return "SR"
    if (dO, dC, dX) == (1, 1, 0):
        return "OX"

    return None


def build_edges(states):
    edges = []

    for a in states:
        for b in states:
            if b["coordinate"] <= a["coordinate"]:
                continue

            transition = classify_transition(a, b)
            if transition is None:
                continue

            # Coordinate convention:
            # OL, C, SR advance by 1;
            # OX increments two counters simultaneously and advances by 2.
            expected_step = 2 if transition == "OX" else 1

            if abs((b["coordinate"] - a["coordinate"]) - expected_step) < 1e-9:
                edges.append((a["state"], b["state"], transition))

    return edges


edges = build_edges(data)

COLOR = {
    "OL": "#2563eb",
    "C":  "#9333ea",
    "SR": "#16a34a",
    "OX": "#dc2626",
}


# ============================================================
# 3. Chemical-potential correction
# ============================================================
# User-defined corrections:
#
# Olation:
#     - Δμ_Mn(OH)2
#
# Condensation:
#     + Δμ_H2O
#
# Oxidation / self-redox:
#     + 1/2 Δμ_H2O
#
# Oxolation:
#     - Δμ_Mn(OH)2 + Δμ_H2O
#
# Because state S(O,C,X) stores cumulative counts:
#
# G = E0 - O*Δμ_Mn(OH)2 + (C + X/2)*Δμ_H2O
# ============================================================
def corrected_energy(s, mu_mn, mu_h2o):
    return (
        s["E0"]
        - s["O"] * mu_mn
        + (s["C"] + 0.5 * s["X"]) * mu_h2o
    )


def edge_arrays(edge_type, mu_mn=0.0, mu_h2o=0.0):
    x, y = [], []

    for a_name, b_name, transition in edges:
        if transition != edge_type:
            continue

        a = smap[a_name]
        b = smap[b_name]

        x += [a["coordinate"], b["coordinate"], None]
        y += [
            corrected_energy(a, mu_mn, mu_h2o),
            corrected_energy(b, mu_mn, mu_h2o),
            None,
        ]

    return x, y


# ============================================================
# 4. Plot
# ============================================================
fig = go.FigureWidget()

for transition in ["OL", "C", "SR", "OX"]:
    x, y = edge_arrays(transition)

    fig.add_scatter(
        x=x,
        y=y,
        mode="lines",
        line=dict(width=2.2, color=COLOR[transition]),
        opacity=0.72,
        name=transition,
        hoverinfo="skip",
    )


# Short horizontal energy levels
level_x = []
level_y = []

marker_x = []
marker_y = []
labels = []
custom = []

for s in data:
    y = corrected_energy(s, 0.0, 0.0)

    level_x += [
        s["coordinate"] - 0.18,
        s["coordinate"] + 0.18,
        None,
    ]
    level_y += [y, y, None]

    marker_x.append(s["coordinate"])
    marker_y.append(y)
    labels.append(s["state"])

    custom.append(
        [
            s["state"],
            s["path"],
            s["E0"],
            s["O"],
            s["C"],
            s["X"],
        ]
    )


fig.add_scatter(
    x=level_x,
    y=level_y,
    mode="lines",
    line=dict(width=3.5, color="#111111"),
    hoverinfo="skip",
    showlegend=False,
)

fig.add_scatter(
    x=marker_x,
    y=marker_y,
    mode="markers+text",
    text=labels,
    textposition="top center",
    marker=dict(size=8, color="#111111"),
    customdata=custom,
    hovertemplate=(
        "<b>%{customdata[0]}</b><br>"
        "Path: %{customdata[1]}<br>"
        "Reaction coordinate: %{x}<br>"
        "Corrected relative energy: %{y:.6f} eV<br>"
        "E0: %{customdata[2]:.6f} eV<br>"
        "O=%{customdata[3]}, "
        "C=%{customdata[4]}, "
        "X=%{customdata[5]}"
        "<extra></extra>"
    ),
    showlegend=False,
)


fig.update_layout(
    title="Chemical-potential dependent free-energy diagram",
    xaxis_title="Reaction coordinate",
    yaxis_title="Relative energy (eV)",
    xaxis=dict(
        dtick=1,
        range=[-0.45, 7.45],
    ),
    legend=dict(
        orientation="h",
        x=0.5,
        xanchor="center",
        y=1.10,
    ),
    template="plotly_white",
    height=780,
)


# ============================================================
# 5. Interactive chemical-potential sliders
# ============================================================
mu_mn = widgets.FloatSlider(
    value=0.0,
    min=-4.0,
    max=0.0,
    step=0.05,
    description="Δμ Mn(OH)2",
    continuous_update=True,
    readout_format=".2f",
    layout=widgets.Layout(width="650px"),
)

mu_h2o = widgets.FloatSlider(
    value=0.0,
    min=-4.0,
    max=0.0,
    step=0.05,
    description="Δμ H2O",
    continuous_update=True,
    readout_format=".2f",
    layout=widgets.Layout(width="650px"),
)


def update(change=None):
    with fig.batch_update():

        # Update reaction-path lines.
        for i, transition in enumerate(["OL", "C", "SR", "OX"]):
            x, y = edge_arrays(
                transition,
                mu_mn.value,
                mu_h2o.value,
            )

            fig.data[i].x = x
            fig.data[i].y = y

        # Update state energy levels and markers.
        level_x = []
        level_y = []
        marker_y = []

        for s in data:
            y = corrected_energy(
                s,
                mu_mn.value,
                mu_h2o.value,
            )

            level_x += [
                s["coordinate"] - 0.18,
                s["coordinate"] + 0.18,
                None,
            ]
            level_y += [y, y, None]
            marker_y.append(y)

        fig.data[4].x = level_x
        fig.data[4].y = level_y
        fig.data[5].y = marker_y

        fig.layout.title = (
            "Chemical-potential dependent free-energy diagram"
            f"<br><sup>"
            f"ΔμMn(OH)2 = {mu_mn.value:.2f} eV, "
            f"ΔμH2O = {mu_h2o.value:.2f} eV"
            f"</sup>"
        )

        fig.layout.yaxis.autorange = True


mu_mn.observe(update, names="value")
mu_h2o.observe(update, names="value")


print(f"Loaded {len(data)} states from: {CSV_FILE}")
print(f"Generated {len(edges)} reaction-network edges.")
display(widgets.VBox([mu_mn, mu_h2o]), fig)
