from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# ============================================================
# Chemical-potential ranges
# 실제 열역학적으로 허용되는 범위로 수정
# ============================================================
MU_MN_MIN = -3.0
MU_MN_MAX = 0.0

MU_H2O_MIN = -2.0
MU_H2O_MAX = 0.0

GRID_SIZE = 121


@st.cache_data
def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    states = pd.read_csv("state_coefficients.csv")
    paths = pd.read_csv("paths.csv")

    required_states = {"state", "C", "a_mn", "a_h2o"}
    required_paths = {"path", "order", "state", "reaction"}

    missing_states = required_states - set(states.columns)
    missing_paths = required_paths - set(paths.columns)

    if missing_states:
        raise ValueError(
            f"state_coefficients.csv에 없는 열: "
            f"{sorted(missing_states)}"
        )

    if missing_paths:
        raise ValueError(
            f"paths.csv에 없는 열: "
            f"{sorted(missing_paths)}"
        )

    if states["state"].duplicated().any():
        duplicated = states.loc[
            states["state"].duplicated(), "state"
        ].tolist()
        raise ValueError(f"중복된 state가 있습니다: {duplicated}")

    return states, paths


def get_route(
    states: pd.DataFrame,
    paths: pd.DataFrame,
    path_name: str,
    mu_mn: float,
    mu_h2o: float,
) -> pd.DataFrame:
    route = (
        paths.loc[paths["path"] == path_name]
        .merge(
            states,
            on="state",
            how="left",
            validate="many_to_one",
        )
        .sort_values("order")
        .reset_index(drop=True)
    )

    if route[["C", "a_mn", "a_h2o"]].isna().any().any():
        missing = route.loc[
            route["C"].isna(), "state"
        ].tolist()
        raise ValueError(
            f"계수 정보가 없는 state가 있습니다: {missing}"
        )

    route["G"] = (
        route["C"]
        + route["a_mn"] * mu_mn
        + route["a_h2o"] * mu_h2o
    )

    # 시작 state를 현재 chemical potential에서 항상 0 eV로 설정
    route["G_rel"] = route["G"] - route.loc[0, "G"]

    # 해당 state로 들어오는 reaction의 ΔG
    route["dG"] = route["G_rel"].diff()

    return route


def get_fixed_y_range(route: pd.DataFrame) -> tuple[float, float]:
    """
    슬라이더를 움직여도 y축이 계속 재조정되지 않도록,
    chemical-potential 영역 전체를 포함하는 y축 범위를 계산.
    선형식이므로 직사각형 영역의 네 모서리만 평가하면 충분함.
    """
    values: list[float] = []

    for mu_mn in (MU_MN_MIN, MU_MN_MAX):
        for mu_h2o in (MU_H2O_MIN, MU_H2O_MAX):
            g = (
                route["C"].to_numpy()
                + route["a_mn"].to_numpy() * mu_mn
                + route["a_h2o"].to_numpy() * mu_h2o
            )

            g_rel = g - g[0]
            values.extend(g_rel.tolist())

    y_min = min(values)
    y_max = max(values)

    span = y_max - y_min
    padding = max(0.30, 0.08 * span)

    return y_min - padding, y_max + padding


def make_energy_diagram(
    route: pd.DataFrame,
    mu_mn: float,
    mu_h2o: float,
    show_all_dg: bool,
) -> go.Figure:
    fig = go.Figure()

    energies = route["G_rel"].to_numpy()
    states = route["state"].tolist()
    reactions = route["reaction"].fillna("").tolist()

    dgs = np.diff(energies)
    max_index = int(np.argmax(dgs))

    half_width = 0.30

    # --------------------------------------------------------
    # State energy levels
    # --------------------------------------------------------
    for i, (state, energy) in enumerate(zip(states, energies)):
        fig.add_trace(
            go.Scatter(
                x=[i - half_width, i + half_width],
                y=[energy, energy],
                mode="lines",
                line={
                    "color": "#222222",
                    "width": 5,
                },
                hovertemplate=(
                    f"<b>{state}</b><br>"
                    f"G = {energy:+.4f} eV"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        )

    # --------------------------------------------------------
    # Reaction connectors
    # --------------------------------------------------------
    for i, dg in enumerate(dgs):
        source = states[i]
        target = states[i + 1]
        reaction = reactions[i + 1]

        is_maximum = i == max_index

        line_color = "#d62728" if is_maximum else "#888888"
        line_width = 4 if is_maximum else 2

        fig.add_trace(
            go.Scatter(
                x=[
                    i + half_width,
                    i + 1 - half_width,
                ],
                y=[
                    energies[i],
                    energies[i + 1],
                ],
                mode="lines",
                line={
                    "color": line_color,
                    "width": line_width,
                },
                hovertemplate=(
                    f"<b>{source} → {target}</b><br>"
                    f"Reaction: {reaction}<br>"
                    f"ΔG = {dg:+.4f} eV"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        )

        if show_all_dg or is_maximum:
            label = (
                f"<b>{reaction}</b><br>"
                f"ΔG = {dg:+.2f} eV"
            )

            fig.add_annotation(
                x=i + 0.5,
                y=(energies[i] + energies[i + 1]) / 2,
                text=label,
                showarrow=False,
                font={
                    "size": 11,
                    "color": line_color,
                },
                bgcolor="rgba(255,255,255,0.8)",
            )

    y_range = get_fixed_y_range(route)

    fig.update_layout(
        title=(
            "Chemical-potential-dependent energy diagram"
            f"<br><sup>"
            f"ΔμMn(OH)₂ = {mu_mn:+.3f} eV, "
            f"ΔμH₂O = {mu_h2o:+.3f} eV"
            f"</sup>"
        ),
        xaxis={
            "title": "Reaction sequence",
            "tickmode": "array",
            "tickvals": list(range(len(states))),
            "ticktext": states,
            "range": [-0.55, len(states) - 0.45],
        },
        yaxis={
            "title": "Relative Gibbs free energy (eV)",
            "range": list(y_range),
            "zeroline": True,
        },
        template="plotly_white",
        hovermode="closest",
        height=600,
        margin={
            "l": 80,
            "r": 30,
            "t": 100,
            "b": 80,
        },
    )

    return fig


def make_dgmax_heatmap(
    route: pd.DataFrame,
    current_mu_mn: float,
    current_mu_h2o: float,
) -> go.Figure:
    mu_mn_values = np.linspace(
        MU_MN_MIN,
        MU_MN_MAX,
        GRID_SIZE,
    )
    mu_h2o_values = np.linspace(
        MU_H2O_MIN,
        MU_H2O_MAX,
        GRID_SIZE,
    )

    mu_mn_grid, mu_h2o_grid = np.meshgrid(
        mu_mn_values,
        mu_h2o_values,
    )

    c = route["C"].to_numpy()[:, None, None]
    a_mn = route["a_mn"].to_numpy()[:, None, None]
    a_h2o = route["a_h2o"].to_numpy()[:, None, None]

    energies = (
        c
        + a_mn * mu_mn_grid[None, :, :]
        + a_h2o * mu_h2o_grid[None, :, :]
    )

    reaction_dg = np.diff(energies, axis=0)

    # 각 chemical-potential 점에서 가장 큰 reaction ΔG
    dg_max = np.max(reaction_dg, axis=0)

    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=mu_mn_values,
            y=mu_h2o_values,
            z=dg_max,
            colorbar={
                "title": "ΔGmax<br>(eV)",
            },
            hovertemplate=(
                "ΔμMn(OH)₂ = %{x:.3f} eV<br>"
                "ΔμH₂O = %{y:.3f} eV<br>"
                "ΔGmax = %{z:.3f} eV"
                "<extra></extra>"
            ),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[current_mu_mn],
            y=[current_mu_h2o],
            mode="markers",
            marker={
                "size": 14,
                "symbol": "x",
                "color": "black",
            },
            name="Current condition",
            hovertemplate=(
                "Current condition<br>"
                f"ΔμMn(OH)₂ = {current_mu_mn:+.3f} eV<br>"
                f"ΔμH₂O = {current_mu_h2o:+.3f} eV"
                "<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Maximum reaction free energy over chemical-potential space",
        xaxis_title="ΔμMn(OH)₂ (eV)",
        yaxis_title="ΔμH₂O (eV)",
        template="plotly_white",
        height=600,
    )

    return fig


def main() -> None:
    st.set_page_config(
        page_title="Mn₃O₄ Energy Landscape",
        layout="wide",
    )

    st.title("Mn₃O₄ chemical-potential-dependent energy landscape")

    try:
        states, paths = load_data()
    except (FileNotFoundError, ValueError) as exc:
        st.error(str(exc))
        st.stop()

    path_names = paths["path"].drop_duplicates().tolist()

    with st.sidebar:
        st.header("Conditions")

        selected_path = st.selectbox(
            "Reaction pathway",
            options=path_names,
        )

        mu_mn = st.slider(
            "ΔμMn(OH)₂ (eV)",
            min_value=MU_MN_MIN,
            max_value=MU_MN_MAX,
            value=(MU_MN_MIN + MU_MN_MAX) / 2,
            step=0.01,
        )

        mu_h2o = st.slider(
            "ΔμH₂O (eV)",
            min_value=MU_H2O_MIN,
            max_value=MU_H2O_MAX,
            value=(MU_H2O_MIN + MU_H2O_MAX) / 2,
            step=0.01,
        )

        show_all_dg = st.checkbox(
            "Show all ΔG labels",
            value=True,
        )

    route = get_route(
        states=states,
        paths=paths,
        path_name=selected_path,
        mu_mn=mu_mn,
        mu_h2o=mu_h2o,
    )

    dgs = route["dG"].dropna().to_numpy()
    max_index = int(np.argmax(dgs))

    source_state = route.loc[max_index, "state"]
    target_state = route.loc[max_index + 1, "state"]
    max_reaction = route.loc[max_index + 1, "reaction"]
    dg_max = dgs[max_index]

    col1, col2, col3 = st.columns(3)

    col1.metric(
        "Maximum step ΔG",
        f"{dg_max:+.3f} eV",
    )

    col2.metric(
        "Bottleneck transition",
        f"{source_state} → {target_state}",
    )

    col3.metric(
        "Reaction type",
        str(max_reaction),
    )

    energy_fig = make_energy_diagram(
        route=route,
        mu_mn=mu_mn,
        mu_h2o=mu_h2o,
        show_all_dg=show_all_dg,
    )

    st.plotly_chart(
        energy_fig,
        width="stretch",
    )

    heatmap_fig = make_dgmax_heatmap(
        route=route,
        current_mu_mn=mu_mn,
        current_mu_h2o=mu_h2o,
    )

    st.plotly_chart(
        heatmap_fig,
        width="stretch",
    )

    display_columns = [
        "order",
        "state",
        "reaction",
        "C",
        "a_mn",
        "a_h2o",
        "G_rel",
        "dG",
    ]

    st.subheader("Current energy table")
    st.dataframe(
        route[display_columns],
        width="stretch",
        hide_index=True,
    )


if __name__ == "__main__":
    main()