#!/usr/bin/env python3
"""Generate capping-agent inputs for combined Mn3O4(011) layer-12/34 states."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LIGAND_ROOT = Path("/home/sowon-desktop/ligand_adsorption")
PROJECT_ROOT = Path("/home/sowon-desktop/ligand-project/EnergeticSpanModel")
DEFAULT_SOURCE_12 = LIGAND_ROOT / "generated_contcar_011_12"
DEFAULT_SOURCE_34 = LIGAND_ROOT / "generated_contcar_011_34"
DEFAULT_OUTCAR = LIGAND_ROOT / "OUTCAR"
DEFAULT_REFERENCE = DEFAULT_SOURCE_34 / "contcar_s322_01.vasp"
DEFAULT_TEMPLATE_INCAR = PROJECT_ROOT / "011_facet_34/34layer/INCAR"
DEFAULT_TEMPLATE_KPOINTS = PROJECT_ROOT / "011_facet_34/34layer/KPOINTS"
DEFAULT_SELECTION = PROJECT_ROOT / "rNets_diagram/plotly/minimum_path_states_011_1234_mu_-4_0.csv"
DEFAULT_OUTPUT = LIGAND_ROOT / "INPUT/011_minimum_path"
SHIFT_34 = (3, 2, 2)


def load_geometry_engine():
    candidates = (
        HERE / "generate_001_ligand_adsorption_inputs.py",
        LIGAND_ROOT / "generate_001_ligand_adsorption_inputs.py",
        PROJECT_ROOT / "rNets_diagram/plotly/generate_001_ligand_adsorption_inputs.py",
    )
    engine_path = next((path for path in candidates if path.is_file()), None)
    if engine_path is None:
        raise FileNotFoundError("Could not locate generate_001_ligand_adsorption_inputs.py")
    spec = importlib.util.spec_from_file_location("ligand_adsorption_geometry_engine", engine_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load geometry engine: {engine_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, engine_path


def parse_last_magnetization_tot(path: Path) -> list[float]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for index, line in enumerate(lines):
        if "magnetization (x)" in line:
            start = index
    if start is None:
        raise ValueError(f"No magnetization (x) block in {path}")
    row = re.compile(
        r"^\s*(\d+)\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+"
        r"[-+]?\d+\.\d+\s+([-+]?\d+\.\d+)\s*$"
    )
    values = []
    for line in lines[start:]:
        if line.strip().startswith("tot ") and values:
            break
        match = row.match(line)
        if match:
            values.append(float(match.group(2)))
    if not values:
        raise ValueError(f"No ion magnetization rows parsed from {path}")
    return values


def discretize_magmom(value: float) -> float:
    if abs(value) < 0.5:
        return 0.0
    magnitude = 4.4 if abs(value) >= 4.1 else 3.8
    return magnitude if value > 0 else -magnitude


def replace_or_append(text: str, key: str, value: str) -> str:
    replacement = f"{key} = {value}"
    pattern = rf"(?im)^\s*{re.escape(key)}\s*=.*$"
    if re.search(pattern, text):
        return re.sub(pattern, replacement, text)
    return text.rstrip() + "\n" + replacement + "\n"


def write_incar(
    engine,
    template: Path,
    destination: Path,
    slab,
    ligand,
    mn_magmom: list[float],
    final_symbols: list[str],
) -> None:
    text = template.read_text(encoding="utf-8", errors="replace")
    if len(mn_magmom) != slab.symbols.count("Mn"):
        raise ValueError("Mn MAGMOM count does not match slab Mn count")
    mn_cursor = 0
    records: list[tuple[str, float]] = []
    for symbol in slab.symbols:
        if symbol == "Mn":
            records.append((symbol, mn_magmom[mn_cursor]))
            mn_cursor += 1
        else:
            records.append((symbol, 0.0))
    records.extend((symbol, 0.0) for symbol in ligand.symbols)
    ordered = [moment for species in engine.SPECIES_ORDER for symbol, moment in records if symbol == species]
    if len(ordered) != len(final_symbols):
        raise ValueError("Final MAGMOM construction failed")
    text = replace_or_append(text, "MAGMOM", engine.compress_magmom(ordered))
    text = replace_or_append(text, "NCORE", "8")
    species = [symbol for symbol in engine.SPECIES_ORDER if symbol in final_symbols]
    text = replace_or_append(text, "LDAUL", " ".join("2" if s == "Mn" else "-1" for s in species))
    text = replace_or_append(text, "LDAUU", " ".join("3.900" if s == "Mn" else "0.000" for s in species))
    text = replace_or_append(text, "LDAUJ", " ".join("0.000" for _ in species))
    destination.write_text(text.rstrip() + "\n", encoding="utf-8")


def parse_source_name(path: Path) -> tuple[tuple[int, int, int], int]:
    match = re.fullmatch(r"contcar_s(\d)(\d)(\d)_(\d+)\.vasp", path.name, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Unexpected source filename: {path.name}")
    x, y, z, variant = map(int, match.groups())
    return (x, y, z), variant


def state_label(values: tuple[int, int, int]) -> str:
    return "S" + "".join(str(value) for value in values)


def collect_sources(source_12: Path, source_34: Path) -> list[dict]:
    records = []
    for source_set, root, shift in (
        ("011_12", source_12, (0, 0, 0)),
        ("011_34", source_34, SHIFT_34),
    ):
        for path in sorted(root.glob("contcar_s*.vasp")):
            local, variant = parse_source_name(path)
            combined = tuple(value + delta for value, delta in zip(local, shift))
            records.append(
                {
                    "source_set": source_set,
                    "source_path": path,
                    "local_state": state_label(local),
                    "combined_state": state_label(combined),
                    "variant": f"variant_{variant:02d}",
                    "state_shift": shift,
                }
            )
    labels = [record["combined_state"] for record in records]
    if len(labels) != len(set(labels)):
        raise ValueError("Combined state labels are not unique")
    return records


def optimize_chelate_flexible(
    engine,
    ligand,
    oxygen_indices: tuple[int, int],
    carbon_index: int,
    slab,
    mn_index: int,
    midpoint_distance: float,
) -> tuple[np.ndarray, float, float, float, float]:
    try:
        positions, azimuth, distance = engine.optimize_chelate(
            ligand, oxygen_indices, carbon_index, slab, mn_index, midpoint_distance
        )
        return positions, azimuth, distance, 0.0, 0.0
    except ValueError:
        pass
    normal = engine.normal_from_cell(slab.cell)
    lateral_x = slab.cell[0] - np.dot(slab.cell[0], normal) * normal
    lateral_x /= np.linalg.norm(lateral_x)
    lateral_y = np.cross(normal, lateral_x)
    for trial_distance in np.arange(midpoint_distance, midpoint_distance + 0.81, 0.10):
        candidates = []
        for azimuth in range(0, 360, 15):
            base = engine.place_chelate(
                ligand, oxygen_indices, carbon_index, slab, mn_index,
                float(trial_distance), float(azimuth),
            )
            for anchor_tilt in (10.0, 20.0, 30.0, 40.0):
                theta = math.radians(anchor_tilt)
                for anchor_azimuth in range(0, 360, 30):
                    phi = math.radians(anchor_azimuth)
                    shift = trial_distance * (
                        (math.cos(theta) - 1.0) * normal
                        + math.sin(theta)
                        * (math.cos(phi) * lateral_x + math.sin(phi) * lateral_y)
                    )
                    positions = base + shift
                    clearance = min(
                        engine.minimum_ligand_slab_distance(positions, slab),
                        engine.minimum_ligand_periodic_nonbonded_distance(
                            positions, ligand, slab.cell
                        ),
                    )
                    candidates.append(
                        (
                            clearance, positions, float(azimuth),
                            float(trial_distance), anchor_tilt, float(anchor_azimuth),
                        )
                    )
        best = max(candidates, key=lambda item: item[0])
        if best[0] >= engine.MIN_INTERATOMIC_DISTANCE:
            return best[1], best[2], round(best[3], 4), best[4], best[5]
    raise ValueError(f"No clash-free flexible OO-chelate geometry at Mn {mn_index + 1}")


def build_case(
    engine,
    slab,
    ligand,
    positions: np.ndarray,
    source: dict,
    case_dir: Path,
    metadata: dict,
    mn_magmom: list[float],
    template_incar: Path,
    template_kpoints: Path,
    potcar: Path,
    outcar: Path,
) -> dict:
    slab_clearance = engine.minimum_ligand_slab_distance(positions, slab)
    image_clearance = engine.minimum_ligand_periodic_nonbonded_distance(positions, ligand, slab.cell)
    if min(slab_clearance, image_clearance) < engine.MIN_INTERATOMIC_DISTANCE:
        raise ValueError(
            f"Final geometry clash: slab={slab_clearance:.3f}, image={image_clearance:.3f} A"
        )
    combined = engine.Poscar(
        comment=(
            f"Mn3O4(011) {source['combined_state']} [{source['source_set']} "
            f"{source['local_state']}] + {metadata['ligand']} {metadata['configuration']}"
        ),
        cell=slab.cell.copy(),
        symbols=slab.symbols + ligand.symbols,
        positions=np.vstack((slab.positions, positions)),
        flags=slab.flags + [("T", "T", "T")] * len(ligand.symbols),
    )
    case_dir.mkdir(parents=True, exist_ok=False)
    engine.write_poscar(case_dir / "POSCAR", combined)
    write_incar(engine, template_incar, case_dir / "INCAR", slab, ligand, mn_magmom, combined.symbols)
    shutil.copy2(template_kpoints, case_dir / "KPOINTS")
    species = [symbol for symbol in engine.SPECIES_ORDER if symbol in combined.symbols]
    expected_potcar = ["Mn_pv" if symbol == "Mn" else symbol for symbol in species]
    actual_potcar = engine.potcar_species(potcar)
    if actual_potcar != expected_potcar:
        raise ValueError(f"POTCAR mismatch: expected {expected_potcar}, found {actual_potcar}")
    shutil.copy2(potcar, case_dir / "POTCAR")
    (case_dir / "POTCAR.spec").write_text("\n".join(expected_potcar) + "\n", encoding="utf-8")
    result = {
        **metadata,
        "state": source["combined_state"],
        "source_set": source["source_set"],
        "source_local_state": source["local_state"],
        "source_variant": source["variant"],
        "source_state_shift": list(source["state_shift"]),
        "source_structure": str(source["source_path"]),
        "minimum_path_selection": source["minimum_path_selection"],
        "cooptimal_grid_points": source["cooptimal_grid_points"],
        "representative_grid_points": source["representative_grid_points"],
        "magmom_reference_outcar": str(outcar),
        "magmom_mapping": "Mn sequence index -> same Mn ion index in reference OUTCAR",
        "slab_mn_magmom": mn_magmom,
        "slab_atoms": len(slab.symbols),
        "ligand_atoms": len(ligand.symbols),
        "total_atoms": len(combined.symbols),
        "species_order": species,
        "ligand_magmom": 0.0,
        "minimum_ligand_slab_distance_A": round(slab_clearance, 4),
        "minimum_ligand_periodic_nonbonded_distance_A": round(image_clearance, 4),
        "potcar_source": str(potcar),
    }
    (case_dir / "adsorption_metadata.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-12", type=Path, default=DEFAULT_SOURCE_12)
    parser.add_argument("--source-34", type=Path, default=DEFAULT_SOURCE_34)
    parser.add_argument("--outcar", type=Path, default=DEFAULT_OUTCAR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--template-incar", type=Path, default=DEFAULT_TEMPLATE_INCAR)
    parser.add_argument("--template-kpoints", type=Path, default=DEFAULT_TEMPLATE_KPOINTS)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--ligand-root", type=Path, default=LIGAND_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Output exists; refusing to overwrite: {args.output}")
    required = (
        args.source_12, args.source_34, args.outcar, args.reference,
        args.template_incar, args.template_kpoints,
        args.selection,
        args.ligand_root / "C5-amine.vasp", args.ligand_root / "C5-carboxyl.vasp",
        args.ligand_root / "POTCAR_amine", args.ligand_root / "POTCAR_carboxyl",
    )
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    engine, engine_path = load_geometry_engine()
    reference = engine.read_poscar(args.reference)
    raw_magnetization = parse_last_magnetization_tot(args.outcar)
    if len(raw_magnetization) != len(reference.symbols):
        raise ValueError(
            f"OUTCAR magnetization rows ({len(raw_magnetization)}) != reference atoms "
            f"({len(reference.symbols)})"
        )
    reference_mn_raw = [value for value, symbol in zip(raw_magnetization, reference.symbols) if symbol == "Mn"]
    reference_mn_magmom = [discretize_magmom(value) for value in reference_mn_raw]
    all_sources = collect_sources(args.source_12, args.source_34)
    source_by_state = {source["combined_state"]: source for source in all_sources}
    with args.selection.open(newline="", encoding="utf-8-sig") as handle:
        selection_rows = list(csv.DictReader(handle))
    sources = []
    for row in selection_rows:
        state = row["state"]
        if state not in source_by_state:
            raise ValueError(f"Selected state has no source structure: {state}")
        source = dict(source_by_state[state])
        source["minimum_path_selection"] = row["selection"]
        source["cooptimal_grid_points"] = int(row["cooptimal_grid_points"])
        source["representative_grid_points"] = int(row["representative_grid_points"])
        sources.append(source)
    if len(sources) != 36:
        raise ValueError(f"Expected 36 selected source structures, found {len(sources)}")

    amine = engine.read_poscar(args.ligand_root / "C5-amine.vasp")
    carboxyl = engine.read_poscar(args.ligand_root / "C5-carboxyl.vasp")
    amine_n = [index for index, symbol in enumerate(amine.symbols) if symbol == "N"]
    carboxyl_o = [index for index, symbol in enumerate(carboxyl.symbols) if symbol == "O"]
    carboxyl_c = [index for index, symbol in enumerate(carboxyl.symbols) if symbol == "C"]
    if len(amine_n) != 1 or len(carboxyl_o) != 2:
        raise ValueError("Unexpected ligand anchor atom counts")
    o_midpoint = 0.5 * (carboxyl.positions[carboxyl_o[0]] + carboxyl.positions[carboxyl_o[1]])
    carboxyl_center_c = min(
        carboxyl_c, key=lambda index: np.linalg.norm(carboxyl.positions[index] - o_midpoint)
    )

    args.output.mkdir(parents=True)
    manifest = []
    for source in sources:
        slab = engine.read_poscar(source["source_path"])
        mn_count = slab.symbols.count("Mn")
        if mn_count > len(reference_mn_magmom):
            raise ValueError(f"{source['source_path']} has more Mn than the OUTCAR reference")
        mn_magmom = reference_mn_magmom[:mn_count]
        state = source["combined_state"]

        amine_geometries = []
        amine_longer_distance = []
        for site in engine.choose_surface_mn(slab, engine.BIND_DISTANCE_N_MN):
            base_angle = (0.0, 120.0, 240.0)[len(amine_geometries)]
            try:
                geometry = engine.optimize_monodentate(
                    amine, amine_n[0], slab, site, engine.BIND_DISTANCE_N_MN, base_angle
                )
            except ValueError:
                continue
            if geometry[2] <= engine.BIND_DISTANCE_N_MN + 1e-6:
                amine_geometries.append((site, geometry))
            else:
                amine_longer_distance.append((site, geometry))
            if len(amine_geometries) == 3:
                break
        if len(amine_geometries) < 3:
            amine_longer_distance.sort(key=lambda item: item[1][2])
            amine_geometries.extend(
                amine_longer_distance[: 3 - len(amine_geometries)]
            )
        if len(amine_geometries) < 3:
            raise ValueError(f"{state}: only {len(amine_geometries)} distinct Mn-N sites found")
        for number, (site, geometry) in enumerate(amine_geometries, start=1):
            positions, roll, distance, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth = geometry
            config = f"N_site{number:02d}"
            manifest.append(
                build_case(
                    engine, slab, amine, positions, source,
                    args.output / state / "amine" / config,
                    {
                        "ligand": "C5-amine", "configuration": config,
                        "binding_mode": "N_monodentate", "target_mn_indices_1based": [site + 1],
                        "nominal_binding_distance_A": engine.BIND_DISTANCE_N_MN,
                        "actual_binding_distance_A": distance, "roll_deg": roll,
                        "backbone_tilt_deg": tilt, "backbone_tilt_azimuth_deg": tilt_azimuth,
                        "anchor_tilt_deg": anchor_tilt, "anchor_azimuth_deg": anchor_azimuth,
                    },
                    mn_magmom, args.template_incar, args.template_kpoints,
                    args.ligand_root / "POTCAR_amine", args.outcar,
                )
            )

        carboxyl_geometries = []
        carboxyl_sites = engine.choose_surface_mn(slab, engine.BIND_DISTANCE_O_MN)
        for site in carboxyl_sites:
            number = len(carboxyl_geometries)
            anchor = carboxyl_o[number % 2]
            base_angle = (30.0, 150.0)[number]
            try:
                geometry = engine.optimize_monodentate(
                    carboxyl, anchor, slab, site, engine.BIND_DISTANCE_O_MN, base_angle
                )
            except ValueError:
                continue
            carboxyl_geometries.append((site, anchor, geometry))
            if len(carboxyl_geometries) == 2:
                break
        if len(carboxyl_geometries) < 2:
            raise ValueError(f"{state}: only {len(carboxyl_geometries)} distinct Mn-O sites found")
        for number, (site, anchor, geometry) in enumerate(carboxyl_geometries, start=1):
            positions, roll, distance, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth = geometry
            config = f"O_site{number:02d}"
            manifest.append(
                build_case(
                    engine, slab, carboxyl, positions, source,
                    args.output / state / "carboxyl" / config,
                    {
                        "ligand": "C5-carboxyl", "configuration": config,
                        "binding_mode": "O_monodentate", "target_mn_indices_1based": [site + 1],
                        "anchor_ligand_index_1based": anchor + 1,
                        "nominal_binding_distance_A": engine.BIND_DISTANCE_O_MN,
                        "actual_binding_distance_A": distance, "roll_deg": roll,
                        "backbone_tilt_deg": tilt, "backbone_tilt_azimuth_deg": tilt_azimuth,
                        "anchor_tilt_deg": anchor_tilt, "anchor_azimuth_deg": anchor_azimuth,
                    },
                    mn_magmom, args.template_incar, args.template_kpoints,
                    args.ligand_root / "POTCAR_carboxyl", args.outcar,
                )
            )

        mono_sites = {site for site, _, _ in carboxyl_geometries}
        chelate_candidates = [site for site in carboxyl_sites if site not in mono_sites]
        chelate = None
        for site in chelate_candidates:
            try:
                geometry = optimize_chelate_flexible(
                    engine,
                    carboxyl, (carboxyl_o[0], carboxyl_o[1]), carboxyl_center_c,
                    slab, site, engine.BIND_DISTANCE_O_MN,
                )
            except ValueError:
                continue
            positions = geometry[0]
            o_mn_distances = [
                engine.distance_pbc_xy(positions[index], slab.positions[site], slab.cell)
                for index in carboxyl_o
            ]
            if not all(1.8 <= distance <= 2.7 for distance in o_mn_distances):
                continue
            chelate = (site, geometry, o_mn_distances)
            break
        if chelate is not None:
            site, geometry, o_mn_distances = chelate
            positions, azimuth, midpoint_distance, anchor_tilt, anchor_azimuth = geometry
            config = "OO_chelate01"
            manifest.append(
                build_case(
                    engine, slab, carboxyl, positions, source,
                    args.output / state / "carboxyl" / config,
                    {
                        "ligand": "C5-carboxyl", "configuration": config,
                        "binding_mode": "OO_chelate", "target_mn_indices_1based": [site + 1],
                        "nominal_midpoint_distance_A": engine.BIND_DISTANCE_O_MN,
                        "actual_midpoint_distance_A": midpoint_distance,
                        "actual_O_Mn_distances_A": [round(value, 4) for value in o_mn_distances],
                        "azimuth_deg": azimuth, "anchor_tilt_deg": anchor_tilt,
                        "anchor_azimuth_deg": anchor_azimuth,
                    },
                    mn_magmom, args.template_incar, args.template_kpoints,
                    args.ligand_root / "POTCAR_carboxyl", args.outcar,
                )
            )
        else:
            third = None
            for site in chelate_candidates:
                for anchor in carboxyl_o:
                    try:
                        geometry = engine.optimize_monodentate(
                            carboxyl, anchor, slab, site,
                            engine.BIND_DISTANCE_O_MN, 270.0,
                        )
                    except ValueError:
                        continue
                    third = (site, anchor, geometry)
                    break
                if third is not None:
                    break
            if third is None:
                raise ValueError(f"{state}: no distinct third carboxyl adsorption site found")
            site, anchor, geometry = third
            positions, roll, distance, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth = geometry
            config = "O_site03"
            manifest.append(
                build_case(
                    engine, slab, carboxyl, positions, source,
                    args.output / state / "carboxyl" / config,
                    {
                        "ligand": "C5-carboxyl", "configuration": config,
                        "binding_mode": "O_monodentate", "target_mn_indices_1based": [site + 1],
                        "anchor_ligand_index_1based": anchor + 1,
                        "nominal_binding_distance_A": engine.BIND_DISTANCE_O_MN,
                        "actual_binding_distance_A": distance, "roll_deg": roll,
                        "backbone_tilt_deg": tilt, "backbone_tilt_azimuth_deg": tilt_azimuth,
                        "anchor_tilt_deg": anchor_tilt, "anchor_azimuth_deg": anchor_azimuth,
                        "fallback_reason": "No distinct OO-chelate geometry with both O-Mn distances in 1.8-2.7 A",
                    },
                    mn_magmom, args.template_incar, args.template_kpoints,
                    args.ligand_root / "POTCAR_carboxyl", args.outcar,
                )
            )

    columns = (
        "state", "source_set", "source_local_state", "source_variant", "source_state_shift",
        "source_structure", "minimum_path_selection", "cooptimal_grid_points",
        "representative_grid_points", "ligand", "configuration", "binding_mode",
        "target_mn_indices_1based", "slab_atoms", "ligand_atoms", "total_atoms",
        "species_order", "minimum_ligand_slab_distance_A",
        "minimum_ligand_periodic_nonbonded_distance_A",
    )
    with (args.output / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for item in manifest:
            row = dict(item)
            for key in ("source_state_shift", "target_mn_indices_1based", "species_order"):
                row[key] = " ".join(map(str, row[key]))
            writer.writerow(row)

    oo_count = sum(item["binding_mode"] == "OO_chelate" for item in manifest)
    fallback_states = sorted(
        item["state"] for item in manifest if item["configuration"] == "O_site03"
    )
    selected_states_text = ",".join(source["combined_state"] for source in sources)
    readme = f"""# Mn3O4(011) combined ligand adsorption inputs

Generated cases: {len(manifest)} ({len(sources)} selected states x 6 configurations)

- Layer-12 source states retain their original state labels.
- Layer-34 source labels are shifted by +(3,2,2), e.g. local S000 -> combined S322.
- Flat hierarchy: `<combined-state>/<ligand>/<configuration>`.
- Each state has three amine cases and three carboxyl cases.
- Carboxyl uses two O-monodentate cases plus an OO-chelate when both O-Mn
  distances are 1.8-2.7 A; otherwise a third distinct O-monodentate site is used.
- OO-chelate cases: {oo_count}; O_site03 fallbacks: {len(fallback_states)}
  ({','.join(fallback_states) if fallback_states else 'none'}).
- Mn MAGMOM values come from the last magnetization block of `{args.outcar}`,
  discretized to +/-3.8 or +/-4.4 and mapped by Mn sequence index.
- Slab O/H and every ligand atom have MAGMOM 0.0.
- Geometry engine: `{engine_path}`.
- Source variant and local state remain recorded in metadata and manifest.
- State selection: `{args.selection}`; union of minimum energetic-span paths on
  the 81 x 81 Mn(OH)2/H2O chemical-potential grid from -4 to 0 eV.
- Selected states: {selected_states_text}
"""
    (args.output / "README.md").write_text(readme, encoding="utf-8")
    print(f"Generated {len(manifest)} cases in {args.output}")


if __name__ == "__main__":
    main()
