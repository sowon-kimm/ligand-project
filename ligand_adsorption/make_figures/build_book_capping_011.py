#!/usr/bin/env python3
"""Build the [011] pure/carboxyl/amine relative-energy table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_PROJECT = Path("/home/sowon-desktop/ligand-project")
DEFAULT_WORKBOOK = DEFAULT_PROJECT / "ligand_adsorption/make_figures/ligand_011.xlsx"
DEFAULT_PURE = (
    DEFAULT_PROJECT
    / "EnergeticSpanModel/rNets_diagram/plotly/csv_outputs"
    / "minimum_path_states_011_1234_mu_-4_0.csv"
)
DEFAULT_OUTPUT = DEFAULT_PROJECT / "ligand_adsorption/make_figures/book_capping_011.csv"

# Same DFT references and reservoir corrections used for book_capping.CSV (001).
PURE_S000_E0 = -266.850780
LIGAND_E0 = {"carboxyl": -90.650006, "amine": -101.803290}
RESERVOIR_E0 = {"O": 30.405611, "C": -14.881036, "X": -4.970412}


def parse_state(state: str) -> tuple[int, int, int]:
    if re.fullmatch(r"S\d{3}", state) is None:
        raise ValueError(f"Invalid state label: {state!r}")
    return tuple(int(value) for value in state[1:])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand-workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--pure-energy", type=Path, default=DEFAULT_PURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--pure-s000-e0", type=float, default=PURE_S000_E0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ligand = pd.read_excel(args.ligand_workbook, sheet_name="lowest_energy").iloc[:, :5]
    ligand = ligand[["ligand", "state", "e0"]].copy()
    ligand["state"] = ligand["state"].astype(str).str.strip()
    ligand["ligand"] = ligand["ligand"].astype(str).str.strip().str.lower()
    ligand["e0"] = pd.to_numeric(ligand["e0"], errors="raise")

    duplicate = ligand.duplicated(["ligand", "state"], keep=False)
    if duplicate.any():
        rows = ligand.loc[duplicate, ["ligand", "state"]].to_dict("records")
        raise ValueError(f"Duplicate ligand/state minima: {rows}")

    capped_e0 = ligand.pivot(index="state", columns="ligand", values="e0")
    pure = pd.read_csv(args.pure_energy)
    pure["state"] = pure["state"].astype(str).str.strip()
    pure["relative_E"] = pd.to_numeric(pure["relative_E"], errors="raise")

    required_ligands = list(LIGAND_E0)
    missing_columns = set(required_ligands) - set(capped_e0.columns)
    if missing_columns:
        raise ValueError(f"Missing ligand columns: {sorted(missing_columns)}")

    joined = pure[["state", "relative_E"]].set_index("state").join(capped_e0)
    incomplete = joined[joined[required_ligands].isna().any(axis=1)]
    complete = joined.dropna(subset=required_ligands).reset_index()

    counts = complete["state"].map(parse_state)
    complete[["O", "C", "X"]] = pd.DataFrame(counts.tolist(), index=complete.index)
    complete["pure_x"] = complete[["O", "C", "X"]].sum(axis=1)
    complete["capping_x"] = complete["pure_x"] + 0.5
    reservoir = sum(complete[key] * value for key, value in RESERVOIR_E0.items())

    for ligand_name, ligand_reference in LIGAND_E0.items():
        complete[f"{ligand_name}_relative"] = (
            complete[ligand_name]
            - args.pure_s000_e0
            + reservoir
            - ligand_reference
        )

    output = complete[
        [
            "state",
            "pure_x",
            "relative_E",
            "capping_x",
            "carboxyl_relative",
            "amine_relative",
        ]
    ].rename(
        columns={
            "relative_E": "pure",
            "carboxyl_relative": "carboxyl",
            "amine_relative": "amine",
        }
    )
    output = output.sort_values(["pure_x", "state"]).reset_index(drop=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, float_format="%.6f")

    print(f"Saved {len(output)} complete states to: {args.output}")
    if not incomplete.empty:
        details = []
        for state, row in incomplete.iterrows():
            missing = [name for name in required_ligands if pd.isna(row[name])]
            details.append(f"{state} ({'/'.join(missing)})")
        print("Excluded states without complete ligand energies: " + ", ".join(details))


if __name__ == "__main__":
    main()
