import json
from pathlib import Path


TEMPLATE_FILE = Path(__file__).with_name("standalone_plotly_controls.js")


def _write_slider_html(fig, output_file, payload):
    script = TEMPLATE_FILE.read_text(encoding="utf-8").replace(
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
        post_script=script,
    )


def write_network_html(fig, output_file, states, edges, barrier_threshold):
    _write_slider_html(
        fig,
        output_file,
        {
            "mode": "network",
            "states": states,
            "edges": edges,
            "barrierThreshold": barrier_threshold,
            "transitions": ["OL", "C", "SR", "OX"],
            "levelHalfWidth": 0.18,
        },
    )


def write_best_path_html(
    fig,
    output_file,
    states,
    all_paths,
    tie_tolerance,
    level_half_width,
    show_all_cooptimal,
):
    _write_slider_html(
        fig,
        output_file,
        {
            "mode": "bestPath",
            "states": states,
            "allPaths": all_paths,
            "tieTolerance": tie_tolerance,
            "transitions": ["OL", "C", "SR", "OX"],
            "levelHalfWidth": level_half_width,
            "showAllCooptimal": show_all_cooptimal,
        },
    )
