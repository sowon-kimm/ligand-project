from __future__ import annotations

import csv
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
from ase.constraints import FixAtoms
from ase.io import read, write
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE_DIR = HERE
REFERENCE_POSCAR = SOURCE_DIR / "contcar_s000_01.vasp"
OUTCAR_MAGNETIZATION = SOURCE_DIR / "OUTCAR_322"
TEMPLATE_ROOT = ROOT / "001_facet" / "generated_state_inputs_001"
OUTPUT_ROOT = ROOT / "generated_001_contcar_initial_inputs"

BOTTOM_Z_CUTOFF = 14.0
SLAB_MN_COUNT = 12
TEMPLATE_FILES = ("INCAR", "KPOINTS", "POTCAR")

FILENAME_RE = re.compile(r"^contcar_s(?P<state>\d{3})_(?P<variant>\d{2})\.vasp$")


def parse_last_magnetization_tot(outcar_path: Path) -> list[float]:
    lines = outcar_path.read_text(errors="replace").splitlines()
    block_start = None
    for idx, line in enumerate(lines):
        if "magnetization (x)" in line:
            block_start = idx
    if block_start is None:
        raise RuntimeError(f"No magnetization (x) block found in {outcar_path}")

    values: list[float] = []
    row_pattern = re.compile(
        r"^\s*(\d+)\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+([-+]?\d+\.\d+)\s*$"
    )
    for line in lines[block_start:]:
        if line.strip().startswith("tot ") and values:
            break
        match = row_pattern.match(line)
        if match:
            values.append(float(match.group(2)))

    if not values:
        raise RuntimeError(f"No ion magnetization values parsed from {outcar_path}")
    return values


def magmom_from_tot(tot: float) -> float:
    if abs(tot) < 0.5:
        return 0.0
    magnitude = 4.4 if abs(tot) >= 4.1 else 3.8
    return magnitude if tot > 0 else -magnitude


def compress_magmom(values: list[float]) -> str:
    if not values:
        return ""

    parts: list[str] = []
    current = values[0]
    count = 1

    def fmt(value: float) -> str:
        return f"{value:.4f}"

    for value in values[1:]:
        if value == current:
            count += 1
            continue
        parts.append(f"{count}*{fmt(current)}" if count > 1 else fmt(current))
        current = value
        count = 1

    parts.append(f"{count}*{fmt(current)}" if count > 1 else fmt(current))
    return " ".join(parts)


def parse_contcar_name(path: Path) -> tuple[str, int]:
    match = FILENAME_RE.match(path.name)
    if not match:
        raise RuntimeError(f"Unexpected CONTCAR filename: {path.name}")
    return f"S{match.group('state')}", int(match.group("variant"))


def atom_ranks_by_symbol(symbols: list[str]) -> list[int]:
    counts: defaultdict[str, int] = defaultdict(int)
    ranks: list[int] = []
    for symbol in symbols:
        ranks.append(counts[symbol])
        counts[symbol] += 1
    return ranks


def reference_scaled_positions_by_symbol_rank():
    reference = read(REFERENCE_POSCAR)
    symbols = reference.get_chemical_symbols()
    ranks = atom_ranks_by_symbol(symbols)
    scaled = reference.get_scaled_positions(wrap=False)

    by_rank = {}
    bottom_original = set()
    for idx, (symbol, rank) in enumerate(zip(symbols, ranks)):
        key = (symbol, rank)
        by_rank[key] = scaled[idx].copy()
        if reference.positions[idx, 2] < BOTTOM_Z_CUTOFF:
            bottom_original.add(key)

    return reference, by_rank, bottom_original


def apply_bottom_fix(atoms, reference_scaled_by_rank, bottom_original_keys):
    symbols = atoms.get_chemical_symbols()
    ranks = atom_ranks_by_symbol(symbols)
    scaled = atoms.get_scaled_positions(wrap=False)
    fixed_indices: list[int] = []

    for idx, (symbol, rank) in enumerate(zip(symbols, ranks)):
        key = (symbol, rank)
        if key in bottom_original_keys:
            scaled[idx] = reference_scaled_by_rank[key]
            fixed_indices.append(idx)
        elif atoms.positions[idx, 2] < BOTTOM_Z_CUTOFF:
            fixed_indices.append(idx)

    atoms.set_scaled_positions(scaled)
    atoms.set_constraint(FixAtoms(indices=fixed_indices))
    return fixed_indices


def load_template_metadata(template_dir: Path) -> dict:
    metadata_path = template_dir / "state_metadata.json"
    if not metadata_path.exists():
        return {}
    return json.loads(metadata_path.read_text())


def build_mn_magmom_maps() -> tuple[list[float], dict[int, float]]:
    outcar_values = parse_last_magnetization_tot(OUTCAR_MAGNETIZATION)
    s322_atoms = read(SOURCE_DIR / "contcar_s322_02.vasp")
    symbols = s322_atoms.get_chemical_symbols()
    if len(outcar_values) < len(symbols):
        raise RuntimeError(
            f"{OUTCAR_MAGNETIZATION} has {len(outcar_values)} magnetization rows, "
            f"but contcar_s322_02.vasp has {len(symbols)} atoms"
        )

    mn_values = [
        magmom_from_tot(outcar_values[idx])
        for idx, symbol in enumerate(symbols)
        if symbol == "Mn"
    ]
    if len(mn_values) < SLAB_MN_COUNT:
        raise RuntimeError(f"Expected at least {SLAB_MN_COUNT} Mn values from OUTCAR")

    added_values = mn_values[SLAB_MN_COUNT:]
    site_magmom = {site_id: value for site_id, value in enumerate(added_values, start=1)}
    return mn_values[:SLAB_MN_COUNT], site_magmom


def mic_scaled_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    delta = np.array(a) - np.array(b)
    delta[:2] -= np.round(delta[:2])
    return delta


def cart_distance(cell: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(mic_scaled_delta(a, b) @ cell))


def build_mn_references():
    outcar_values = parse_last_magnetization_tot(OUTCAR_MAGNETIZATION)
    s322_atoms = read(SOURCE_DIR / "contcar_s322_02.vasp")
    symbols = s322_atoms.get_chemical_symbols()
    scaled = s322_atoms.get_scaled_positions(wrap=False)
    references = []
    mn_rank = 0

    for atom_idx, symbol in enumerate(symbols):
        if symbol != "Mn":
            continue
        mn_rank += 1
        references.append(
            {
                "atom_index_1based": atom_idx + 1,
                "mn_rank_1based": mn_rank,
                "scaled_position": scaled[atom_idx].copy(),
                "magmom": magmom_from_tot(outcar_values[atom_idx]),
            }
        )

    return s322_atoms.cell.array, references


def magmoms_for_atoms_by_position(atoms, reference_cell: np.ndarray, mn_references: list[dict]):
    symbols = atoms.get_chemical_symbols()
    scaled = atoms.get_scaled_positions(wrap=False)
    mn_indices = [idx for idx, symbol in enumerate(symbols) if symbol == "Mn"]
    values = [0.0] * len(atoms)

    if len(mn_indices) > len(mn_references):
        raise RuntimeError(f"Current structure has {len(mn_indices)} Mn, but only {len(mn_references)} Mn references")

    cost = np.zeros((len(mn_indices), len(mn_references)))
    for row, atom_idx in enumerate(mn_indices):
        for col, reference in enumerate(mn_references):
            cost[row, col] = cart_distance(reference_cell, scaled[atom_idx], reference["scaled_position"])

    rows, cols = linear_sum_assignment(cost)
    matches = []
    for row, col in zip(rows, cols):
        atom_idx = mn_indices[row]
        reference = mn_references[col]
        values[atom_idx] = reference["magmom"]
        matches.append(
            {
                "atom_index_1based": atom_idx + 1,
                "matched_outcar_atom_index_1based": reference["atom_index_1based"],
                "matched_outcar_mn_rank_1based": reference["mn_rank_1based"],
                "distance_angstrom": round(float(cost[row, col]), 6),
                "magmom": reference["magmom"],
            }
        )

    matches.sort(key=lambda item: item["atom_index_1based"])
    return values, matches


def update_incar_magmom(incar_path: Path, magmoms: list[float]) -> None:
    lines = incar_path.read_text().splitlines()
    filtered = []
    for line in lines:
        upper = line.strip().upper()
        if upper.startswith("MAGMOM") or upper.startswith("NCORE"):
            continue
        filtered.append(line.rstrip())

    filtered.append(f"MAGMOM = {compress_magmom(magmoms)}")
    filtered.append("NCORE = 8")
    incar_path.write_text("\n".join(filtered) + "\n")


def copy_templates(template_dir: Path, output_dir: Path) -> None:
    for filename in TEMPLATE_FILES:
        src = template_dir / filename
        if not src.exists():
            raise RuntimeError(f"Missing template file: {src}")
        shutil.copy2(src, output_dir / filename)


def main() -> None:
    _, reference_scaled_by_rank, bottom_original_keys = reference_scaled_positions_by_symbol_rank()
    slab_mn_magmoms, site_magmom = build_mn_magmom_maps()
    reference_cell, mn_references = build_mn_references()

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    for contcar_path in sorted(SOURCE_DIR.glob("contcar_s*.vasp")):
        state, variant = parse_contcar_name(contcar_path)
        template_dir = TEMPLATE_ROOT / state / f"variant_{variant:02d}"
        output_dir = OUTPUT_ROOT / state / f"variant_{variant:02d}"
        output_dir.mkdir(parents=True, exist_ok=True)

        metadata = load_template_metadata(template_dir)
        atoms = read(contcar_path)
        atoms.arrays.pop("momenta", None)
        fixed_indices = apply_bottom_fix(atoms, reference_scaled_by_rank, bottom_original_keys)
        magmoms, magmom_matches = magmoms_for_atoms_by_position(atoms, reference_cell, mn_references)

        copy_templates(template_dir, output_dir)
        write(output_dir / "POSCAR", atoms, format="vasp", vasp5=True, direct=True, sort=False)
        update_incar_magmom(output_dir / "INCAR", magmoms)

        out_metadata = dict(metadata)
        out_metadata.update(
            {
                "source_contcar": str(contcar_path.relative_to(ROOT)),
                "bottom_fix_z_cutoff": BOTTOM_Z_CUTOFF,
                "fixed_atom_count": len(fixed_indices),
                "fixed_atom_indices_1based": [idx + 1 for idx in fixed_indices],
                "magmom_source_outcar": str(OUTCAR_MAGNETIZATION.relative_to(ROOT)),
                "magmom_assignment": "position_match_to_OUTCAR_322_Mn_sites",
                "magmom_matches": magmom_matches,
                "magmom": compress_magmom(magmoms),
            }
        )
        (output_dir / "state_metadata.json").write_text(json.dumps(out_metadata, indent=2) + "\n")

        manifest_rows.append(
            {
                "state": state,
                "variant": variant,
                "atom_count": len(atoms),
                "fixed_atom_count": len(fixed_indices),
                "magmom": compress_magmom(magmoms),
                "source_contcar": str(contcar_path.relative_to(ROOT)),
                "output_dir": str(output_dir.relative_to(ROOT)),
            }
        )

    with (OUTPUT_ROOT / "manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "state",
                "variant",
                "atom_count",
                "fixed_atom_count",
                "magmom",
                "source_contcar",
                "output_dir",
            ],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Wrote {len(manifest_rows)} input sets to {OUTPUT_ROOT}")
    print(f"Slab Mn MAGMOM from OUTCAR: {compress_magmom(slab_mn_magmoms)}")
    print(f"Added Mn site MAGMOM from OUTCAR: {site_magmom}")


if __name__ == "__main__":
    main()
