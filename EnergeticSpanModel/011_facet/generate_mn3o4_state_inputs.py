from __future__ import annotations

import csv
import json
import math
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read, write


HERE = Path(__file__).resolve().parent
START_DIR = HERE / "12layer_magmom" / "0"
FINAL_DIR = HERE / "12layer_3rd_magmom" / "oxooxooxiolaoxi"
OUTPUT_DIR = HERE / "generated_state_inputs"
START_STRUCTURE_FILE = "POSCAR"
FINAL_STRUCTURE_FILE = "POSCAR"

TEMPLATE_FILES = ("INCAR", "KPOINTS", "POTCAR")
TEMPLATE_SOURCE_DIRS: tuple[Path, ...] | None = None
SYMBOL_ORDER = {"Mn": 0, "O": 1, "H": 2}

MAX_O = 3
MAX_C = 2
MAX_X = 2

# From the reference 12-layer ferrimagnetic slab.
# Original slab Mn indices 0-7 are Mn3+, indices 8-11 are Mn2+.
BASE_MN_MAGMOM = {
    0: 3.8,
    1: 3.8,
    2: -3.8,
    3: -3.8,
    4: 3.8,
    5: 3.8,
    6: -3.8,
    7: -3.8,
    8: 4.4,
    9: 4.4,
    10: 4.4,
    11: 4.4,
}


@dataclass
class AtomRecord:
    symbol: str
    scaled_position: np.ndarray
    source: str
    original_index: int | None = None


def mic_scaled_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    delta = np.array(a) - np.array(b)
    delta[:2] -= np.round(delta[:2])
    return delta


def cart_distance(cell: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(mic_scaled_delta(a, b) @ cell))


def match_reference_atoms(start: Atoms, final: Atoms) -> tuple[list[int], list[int]]:
    """Return final Mn/O atom indices that are not matched to the starting slab."""
    from scipy.optimize import linear_sum_assignment

    added_mn: list[int] = []
    added_o: list[int] = []
    cell = final.cell.array

    for symbol, added in (("Mn", added_mn), ("O", added_o)):
        start_indices = [i for i, atom in enumerate(start) if atom.symbol == symbol]
        final_indices = [i for i, atom in enumerate(final) if atom.symbol == symbol]
        start_scaled = start.get_scaled_positions()[start_indices]
        final_scaled = final.get_scaled_positions()[final_indices]

        cost = np.zeros((len(start_indices), len(final_indices)))
        for i, start_pos in enumerate(start_scaled):
            for j, final_pos in enumerate(final_scaled):
                cost[i, j] = cart_distance(cell, final_pos, start_pos)

        _, matched_final_local = linear_sum_assignment(cost)
        matched_final_local = set(matched_final_local)

        for local_idx, final_idx in enumerate(final_indices):
            if local_idx not in matched_final_local:
                added.append(final_idx)

    return added_mn, added_o


def match_h_changes(start: Atoms, final: Atoms, max_slab_h_distance: float = 0.75) -> tuple[list[int], list[int]]:
    """Return removed starting slab H indices and added final-layer H indices."""
    start_h = [i for i, atom in enumerate(start) if atom.symbol == "H"]
    final_h = [i for i, atom in enumerate(final) if atom.symbol == "H"]
    start_scaled = start.get_scaled_positions()
    final_scaled = final.get_scaled_positions()
    cell = final.cell.array

    removed_start_h = []
    for start_idx in start_h:
        nearest = min(
            cart_distance(cell, start_scaled[start_idx], final_scaled[final_idx])
            for final_idx in final_h
        )
        if nearest > max_slab_h_distance:
            removed_start_h.append(start_idx)

    added_final_h = []
    for final_idx in final_h:
        nearest = min(
            cart_distance(cell, final_scaled[final_idx], start_scaled[start_idx])
            for start_idx in start_h
        )
        if nearest > max_slab_h_distance:
            added_final_h.append(final_idx)

    removed_start_h.sort(key=lambda idx: start_scaled[idx, 2], reverse=True)
    added_final_h.sort(key=lambda idx: final_scaled[idx, 2])
    return removed_start_h, added_final_h


def sorted_site_fragments(final: Atoms, added_mn: list[int], added_o: list[int], added_h: list[int]) -> list[dict]:
    """Group the three final-layer Mn sites with nearby final-layer O positions."""
    scaled = final.get_scaled_positions()
    cell = final.cell.array

    # Left-to-right in cartesian x gives stable site labels.
    mn_indices = sorted(added_mn, key=lambda idx: final.positions[idx, 0])

    fragments = []
    for site_id, mn_idx in enumerate(mn_indices, start=1):
        distances = sorted(
            (cart_distance(cell, scaled[o_idx], scaled[mn_idx]), o_idx)
            for o_idx in added_o
        )
        fragments.append(
            {
                "site_id": site_id,
                "mn_index": mn_idx,
                "mn_scaled": scaled[mn_idx].copy(),
                "near_o_indices": [idx for _, idx in distances],
                "near_h_indices": [
                    idx
                    for _, idx in sorted(
                        (cart_distance(cell, scaled[h_idx], scaled[mn_idx]), h_idx)
                        for h_idx in added_h
                    )
                ],
            }
        )

    return fragments


def choose_sites(o_count: int, variant_idx: int, fragments: list[dict]) -> list[dict]:
    if o_count == 0:
        return []

    orders = [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
    ]
    return [fragments[i] for i in orders[variant_idx - 1][:o_count]]


def generated_o_position(site: dict, nth: int, variant_idx: int, cell: np.ndarray) -> np.ndarray:
    angle = 2.0 * math.pi * ((site["site_id"] - 1) / 3.0 + 0.12 * variant_idx + 0.18 * nth)
    cart_offset = np.array(
        [
            1.35 * math.cos(angle),
            1.35 * math.sin(angle),
            1.25 + 0.12 * nth,
        ]
    )
    scaled_offset = cart_offset @ np.linalg.inv(cell)
    position = np.array(site["mn_scaled"]) + scaled_offset
    position[:2] %= 1.0
    return position


def added_o_positions(
    o_required: int,
    selected_sites: list[dict],
    final: Atoms,
    variant_idx: int,
) -> list[np.ndarray]:
    if o_required <= 0:
        return []

    final_scaled = final.get_scaled_positions()
    cell = final.cell.array
    positions: list[np.ndarray] = []
    used_final_o: set[int] = set()

    # Prefer real final-layer O positions close to the selected Mn sites.
    for site in selected_sites:
        for o_idx in site["near_o_indices"]:
            if o_idx in used_final_o:
                continue
            positions.append(final_scaled[o_idx].copy())
            used_final_o.add(o_idx)
            break
        if len(positions) >= o_required:
            return positions[:o_required]

    # Then add the next closest unused final-layer O positions.
    candidates = []
    for site in selected_sites:
        for o_idx in site["near_o_indices"]:
            if o_idx not in used_final_o:
                dist = cart_distance(cell, final_scaled[o_idx], site["mn_scaled"])
                candidates.append((dist, o_idx))
    for _, o_idx in sorted(candidates):
        if o_idx in used_final_o:
            continue
        positions.append(final_scaled[o_idx].copy())
        used_final_o.add(o_idx)
        if len(positions) >= o_required:
            return positions[:o_required]

    # Early olation states can require more O atoms than exist in the final
    # dehydrated layer. Add terminal OH-like O positions as starting guesses.
    nth = 0
    while len(positions) < o_required:
        site = selected_sites[nth % len(selected_sites)]
        positions.append(generated_o_position(site, nth, variant_idx, cell))
        nth += 1

    return positions[:o_required]


def added_h_positions(
    o_positions: list[np.ndarray],
    h_required: int,
    selected_sites: list[dict],
    final: Atoms,
    variant_idx: int,
    cell: np.ndarray,
) -> list[np.ndarray]:
    if h_required <= 0:
        return []
    if not o_positions:
        raise ValueError("Cannot place added H atoms without added O positions.")

    positions = []
    final_scaled = final.get_scaled_positions()
    used_final_h: set[int] = set()

    for site in selected_sites:
        for h_idx in site["near_h_indices"]:
            if h_idx in used_final_h:
                continue
            positions.append(final_scaled[h_idx].copy())
            used_final_h.add(h_idx)
            break
        if len(positions) >= h_required:
            return positions[:h_required]

    candidates = []
    for site in selected_sites:
        for h_idx in site["near_h_indices"]:
            if h_idx not in used_final_h:
                dist = cart_distance(cell, final_scaled[h_idx], site["mn_scaled"])
                candidates.append((dist, h_idx))
    for _, h_idx in sorted(candidates):
        if h_idx in used_final_h:
            continue
        positions.append(final_scaled[h_idx].copy())
        used_final_h.add(h_idx)
        if len(positions) >= h_required:
            return positions[:h_required]

    inv_cell = np.linalg.inv(cell)
    for i in range(h_required - len(positions)):
        o_pos = o_positions[i % len(o_positions)]
        angle = 2.0 * math.pi * ((i + 0.25 * variant_idx) / max(1, h_required))
        cart_offset = np.array(
            [
                0.35 * math.cos(angle),
                0.35 * math.sin(angle),
                0.90,
            ]
        )
        h_pos = np.array(o_pos) + cart_offset @ inv_cell
        h_pos[:2] %= 1.0
        positions.append(h_pos)
    return positions


def top_h_indices(start: Atoms, variant_idx: int, candidates: list[int] | None = None) -> list[int]:
    h_indices = candidates if candidates is not None else [i for i, atom in enumerate(start) if atom.symbol == "H"]
    scaled = start.get_scaled_positions()
    ordered = sorted(h_indices, key=lambda idx: scaled[idx, 2], reverse=True)

    # Rotate the removal order between variants.
    shift = (variant_idx - 1) % len(ordered)
    return ordered[shift:] + ordered[:shift]


def build_state_atoms(
    start: Atoms,
    final: Atoms,
    fragments: list[dict],
    slab_h_removal_order: list[int],
    state: tuple[int, int, int],
    variant_idx: int,
) -> tuple[Atoms, list[float], dict]:
    o_count, c_count, x_count = state
    selected_sites = choose_sites(o_count, variant_idx, fragments)

    added_o_count = 2 * o_count - c_count
    added_h_count = max(0, 2 * o_count - c_count)
    remove_slab_h_count = min(c_count + x_count, len(slab_h_removal_order))

    removed_h = set(top_h_indices(start, variant_idx, slab_h_removal_order)[:remove_slab_h_count])
    start_scaled = start.get_scaled_positions()
    records: list[AtomRecord] = []

    for idx, atom in enumerate(start):
        if idx in removed_h:
            continue
        records.append(
            AtomRecord(
                symbol=atom.symbol,
                scaled_position=start_scaled[idx].copy(),
                source="slab",
                original_index=idx,
            )
        )

    for site in selected_sites:
        records.append(
            AtomRecord(
                symbol="Mn",
                scaled_position=np.array(site["mn_scaled"]).copy(),
                source=f"ads_site_{site['site_id']}",
            )
        )

    o_positions = added_o_positions(added_o_count, selected_sites, final, variant_idx)
    for idx, pos in enumerate(o_positions, start=1):
        records.append(AtomRecord(symbol="O", scaled_position=pos, source=f"added_O_{idx}"))

    h_positions = added_h_positions(
        o_positions=o_positions,
        h_required=added_h_count,
        selected_sites=selected_sites,
        final=final,
        variant_idx=variant_idx,
        cell=start.cell.array,
    )
    for idx, pos in enumerate(h_positions, start=1):
        records.append(AtomRecord(symbol="H", scaled_position=pos, source=f"added_H_{idx}"))

    records.sort(
        key=lambda rec: (
            SYMBOL_ORDER[rec.symbol],
            0 if rec.source == "slab" else 1,
            rec.original_index if rec.original_index is not None else 10_000,
            rec.source,
        )
    )

    atoms = Atoms(
        symbols=[rec.symbol for rec in records],
        scaled_positions=[rec.scaled_position for rec in records],
        cell=start.cell,
        pbc=True,
    )
    atoms.wrap()

    fixed_mask = [pos[2] < 0.17 for pos in atoms.get_scaled_positions()]
    atoms.set_constraint(FixAtoms(mask=fixed_mask))

    magmoms = []
    for rec in records:
        if rec.symbol == "Mn" and rec.source == "slab":
            magmoms.append(BASE_MN_MAGMOM.get(rec.original_index, 0.0))
        else:
            magmoms.append(0.0)

    metadata = {
        "state": f"S{o_count}{c_count}{x_count}",
        "variant": variant_idx,
        "selected_sites": [site["site_id"] for site in selected_sites],
        "added_o_count": added_o_count,
        "added_h_count": added_h_count,
        "removed_slab_h_indices": sorted(removed_h),
        "formula": atoms.get_chemical_formula(),
        "atom_count": len(atoms),
    }
    return atoms, magmoms, metadata


def compress_magmom(values: list[float]) -> str:
    tokens = []
    i = 0
    while i < len(values):
        value = values[i]
        count = 1
        while i + count < len(values) and values[i + count] == value:
            count += 1

        formatted = f"{value:.4f}" if value != 0 else "0.0000"
        if count == 1:
            tokens.append(formatted)
        else:
            tokens.append(f"{count}*{formatted}")
        i += count
    return " ".join(tokens)


def update_incar_magmom(incar_path: Path, magmoms: list[float]) -> None:
    lines = incar_path.read_text().splitlines()
    filtered = []
    for line in lines:
        upper = line.strip().upper()
        if upper.startswith("MAGMOM") or upper.startswith("NCORE"):
            continue
        filtered.append(line)

    filtered.append(f"MAGMOM = {compress_magmom(magmoms)}")
    filtered.append("NCORE = 8")
    incar_path.write_text("\n".join(filtered) + "\n")


def copy_templates(output_dir: Path) -> None:
    template_dirs = TEMPLATE_SOURCE_DIRS or (START_DIR,)
    for filename in TEMPLATE_FILES:
        for template_dir in template_dirs:
            src = template_dir / filename
            if src.exists():
                shutil.copy2(src, output_dir / filename)
                break
        else:
            raise FileNotFoundError(
                f"Could not find template {filename} in "
                + ", ".join(str(path) for path in template_dirs)
            )


def allowed_states() -> list[tuple[int, int, int]]:
    states = []
    for o_count in range(MAX_O + 1):
        for c_count in range(min(MAX_C, o_count) + 1):
            for x_count in range(MAX_X + 1):
                states.append((o_count, c_count, x_count))
    return states


def main() -> None:
    start = read(START_DIR / START_STRUCTURE_FILE)
    final = read(FINAL_DIR / FINAL_STRUCTURE_FILE)
    added_mn, added_o = match_reference_atoms(start, final)
    removed_start_h, added_h = match_h_changes(start, final)
    fragments = sorted_site_fragments(final, added_mn, added_o, added_h)

    OUTPUT_DIR.mkdir(exist_ok=True)
    manifest_rows = []

    for state in allowed_states():
        state_name = f"S{state[0]}{state[1]}{state[2]}"
        for variant_idx in range(1, 4):
            out_dir = OUTPUT_DIR / state_name / f"variant_{variant_idx:02d}"
            out_dir.mkdir(parents=True, exist_ok=True)

            atoms, magmoms, metadata = build_state_atoms(
                start=start,
                final=final,
                fragments=fragments,
                slab_h_removal_order=removed_start_h,
                state=state,
                variant_idx=variant_idx,
            )

            write(out_dir / "POSCAR", atoms, format="vasp", vasp5=True, direct=True, sort=False)
            copy_templates(out_dir)
            update_incar_magmom(out_dir / "INCAR", magmoms)
            (out_dir / "state_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

            manifest_rows.append(
                {
                    "state": state_name,
                    "variant": f"variant_{variant_idx:02d}",
                    "path": str(out_dir.relative_to(HERE)),
                    "formula": metadata["formula"],
                    "atom_count": metadata["atom_count"],
                    "selected_sites": " ".join(map(str, metadata["selected_sites"])),
                    "added_o_count": metadata["added_o_count"],
                    "added_h_count": metadata["added_h_count"],
                    "removed_slab_h_indices": " ".join(map(str, metadata["removed_slab_h_indices"])),
                }
            )

    with (OUTPUT_DIR / "manifest.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"generated states: {len(allowed_states())}")
    print(f"generated input folders: {len(manifest_rows)}")
    print(f"output: {OUTPUT_DIR.relative_to(HERE)}")


if __name__ == "__main__":
    main()
