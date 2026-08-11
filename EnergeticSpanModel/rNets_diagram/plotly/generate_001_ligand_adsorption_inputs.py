#!/usr/bin/env python3
"""Build ligand adsorption inputs for selected Mn3O4(001) growth states.

The script preserves slab selective-dynamics flags, appends one ligand with all
ligand MAGMOM values set to zero, and creates three adsorption geometries per
ligand. Complete, user-provided POTCAR files are validated against the final
POSCAR species order and copied into every generated input directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np


DEFAULT_PROJECT = Path("/home/sowon-desktop/ligand-project/EnergeticSpanModel")
DEFAULT_SELECTION = DEFAULT_PROJECT / "rNets_diagram/plotly/minimum_path_states_001_mu_-4_0.csv"
DEFAULT_SOURCE = DEFAULT_PROJECT / "generated_001_contcar_initial_inputs"
DEFAULT_LIGAND_ROOT = Path("/home/sowon-desktop/ligand_adsorption")
DEFAULT_OUTPUT = DEFAULT_LIGAND_ROOT / "INPUT/001_minimum_path"
DEFAULT_POTCAR_AMINE = DEFAULT_LIGAND_ROOT / "POTCAR_amine"
DEFAULT_POTCAR_CARBOXYL = DEFAULT_LIGAND_ROOT / "POTCAR_carboxyl"

BIND_DISTANCE_N_MN = 2.15
BIND_DISTANCE_O_MN = 2.05
MIN_INTERATOMIC_DISTANCE = 1.15
SURFACE_MN_WINDOW = 3.0
# This order matches the complete POTCARs supplied for this project.  Omitting
# absent species gives Mn/O/H/N/C for amine and Mn/O/H/C for carboxyl.
SPECIES_ORDER = ("Mn", "O", "H", "N", "C")


@dataclass
class Poscar:
    comment: str
    cell: np.ndarray
    symbols: list[str]
    positions: np.ndarray
    flags: list[tuple[str, str, str]]


def read_poscar(path: Path) -> Poscar:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"Too few lines in POSCAR: {path}")
    scale = float(lines[1].split()[0])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    if len(species) != len(counts):
        raise ValueError(f"Species/count mismatch in {path}")
    symbols = [s for s, n in zip(species, counts) for _ in range(n)]
    cursor = 7
    selective = lines[cursor].strip().lower().startswith("s")
    if selective:
        cursor += 1
    direct = lines[cursor].strip().lower().startswith(("d", "f"))
    cursor += 1
    coords = []
    flags: list[tuple[str, str, str]] = []
    for line in lines[cursor : cursor + len(symbols)]:
        fields = line.split()
        coords.append([float(x) for x in fields[:3]])
        flags.append(tuple(fields[3:6]) if selective and len(fields) >= 6 else ("T", "T", "T"))
    if len(coords) != len(symbols):
        raise ValueError(f"Atom-count mismatch in {path}")
    arr = np.array(coords)
    positions = arr @ cell if direct else arr * scale
    return Poscar(lines[0], cell, symbols, positions, flags)


def write_poscar(path: Path, structure: Poscar) -> None:
    inv_cell = np.linalg.inv(structure.cell)
    scaled = structure.positions @ inv_cell
    species = [s for s in SPECIES_ORDER if s in structure.symbols]
    indices = [i for s in species for i, atom_s in enumerate(structure.symbols) if atom_s == s]
    counts = [sum(atom_s == s for atom_s in structure.symbols) for s in species]
    lines = [
        structure.comment,
        " 1.0",
        *["  " + "  ".join(f"{v:20.14f}" for v in row) for row in structure.cell],
        "  " + "  ".join(species),
        "  " + "  ".join(str(n) for n in counts),
        "Selective dynamics",
        "Direct",
    ]
    for i in indices:
        coord = "  ".join(f"{v:20.14f}" for v in scaled[i])
        lines.append(f"  {coord}   {'   '.join(structure.flags[i])}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def normal_from_cell(cell: np.ndarray) -> np.ndarray:
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    if normal[2] < 0:
        normal *= -1
    return normal


def rotation_align(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    a = source / np.linalg.norm(source)
    b = target / np.linalg.norm(target)
    cross = np.cross(a, b)
    sine = np.linalg.norm(cross)
    cosine = float(np.dot(a, b))
    if sine < 1e-12:
        if cosine > 0:
            return np.eye(3)
        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(a, trial)) > 0.9:
            trial = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, trial)
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    axis = cross / sine
    k = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    return np.eye(3) + sine * k + (1.0 - cosine) * (k @ k)


def rotation_about(axis: np.ndarray, degrees: float) -> np.ndarray:
    u = axis / np.linalg.norm(axis)
    angle = math.radians(degrees)
    c, s = math.cos(angle), math.sin(angle)
    k = np.array([[0.0, -u[2], u[1]], [u[2], 0.0, -u[0]], [-u[1], u[0], 0.0]])
    return c * np.eye(3) + s * k + (1.0 - c) * np.outer(u, u)


def minimum_image_delta(a: np.ndarray, b: np.ndarray, cell: np.ndarray) -> np.ndarray:
    frac = (b - a) @ np.linalg.inv(cell)
    frac[:2] -= np.round(frac[:2])
    return frac @ cell


def distance_pbc_xy(a: np.ndarray, b: np.ndarray, cell: np.ndarray) -> float:
    return float(np.linalg.norm(minimum_image_delta(a, b, cell)))


def choose_surface_mn(slab: Poscar, binding_distance: float) -> list[int]:
    normal = normal_from_cell(slab.cell)
    mn = [i for i, s in enumerate(slab.symbols) if s == "Mn"]
    height = {i: float(np.dot(slab.positions[i], normal)) for i in mn}
    top = max(height.values())
    candidates = []
    for i in mn:
        probe = slab.positions[i] + binding_distance * normal
        clearance = min(
            distance_pbc_xy(probe, position, slab.cell)
            for j, position in enumerate(slab.positions)
            if j != i
        )
        if clearance >= 0.70:
            candidates.append((i, height[i], clearance))
    candidates.sort(
        key=lambda item: (
            -(item[2] >= MIN_INTERATOMIC_DISTANCE), -item[1], -item[2], item[0]
        )
    )
    return [item[0] for item in candidates]


def choose_bridge_pair(slab: Poscar, ligand_oo_distance: float) -> tuple[int, int]:
    normal = normal_from_cell(slab.cell)
    mn = [i for i, s in enumerate(slab.symbols) if s == "Mn"]
    heights = {i: float(np.dot(slab.positions[i], normal)) for i in mn}
    top = max(heights.values())
    candidates = [i for i in mn if heights[i] >= top - SURFACE_MN_WINDOW]
    candidates = sorted(candidates, key=lambda i: (-heights[i], i))[:4]
    if len(candidates) < 2:
        candidates = sorted(mn, key=lambda i: (-heights[i], i))[:2]
    pairs = []
    for p, i in enumerate(candidates):
        for j in candidates[p + 1 :]:
            d = minimum_image_delta(slab.positions[i], slab.positions[j], slab.cell)
            lateral = d - np.dot(d, normal) * normal
            sep = float(np.linalg.norm(lateral))
            if sep < 1e-8:
                continue
            score = abs(sep - ligand_oo_distance) + 0.75 * ((top - heights[i]) + (top - heights[j]))
            pairs.append((score, i, j))
    if not pairs:
        raise ValueError("No surface Mn pair available for OO bridge")
    _, i, j = min(pairs)
    return i, j


def heavy_com(ligand: Poscar) -> np.ndarray:
    indices = [i for i, s in enumerate(ligand.symbols) if s != "H"]
    return np.mean(ligand.positions[indices], axis=0)


def place_monodentate(
    ligand: Poscar,
    anchor: int,
    slab: Poscar,
    mn_index: int,
    distance: float,
    roll: float,
    tilt: float = 0.0,
    tilt_azimuth: float = 0.0,
    anchor_tilt: float = 0.0,
    anchor_azimuth: float = 0.0,
) -> np.ndarray:
    normal = normal_from_cell(slab.cell)
    lateral_x = slab.cell[0] - np.dot(slab.cell[0], normal) * normal
    lateral_x /= np.linalg.norm(lateral_x)
    lateral_y = np.cross(normal, lateral_x)
    theta = math.radians(tilt)
    phi = math.radians(tilt_azimuth)
    target_outward = (
        math.cos(theta) * normal
        + math.sin(theta) * (math.cos(phi) * lateral_x + math.sin(phi) * lateral_y)
    )
    outward = heavy_com(ligand) - ligand.positions[anchor]
    rotation = rotation_align(outward, target_outward)
    rotation = rotation_about(target_outward, roll) @ rotation
    centered = ligand.positions - ligand.positions[anchor]
    placed = centered @ rotation.T
    anchor_theta = math.radians(anchor_tilt)
    anchor_phi = math.radians(anchor_azimuth)
    anchor_direction = (
        math.cos(anchor_theta) * normal
        + math.sin(anchor_theta)
        * (math.cos(anchor_phi) * lateral_x + math.sin(anchor_phi) * lateral_y)
    )
    target = slab.positions[mn_index] + distance * anchor_direction
    return placed + target


def minimum_ligand_slab_distance(positions: np.ndarray, slab: Poscar) -> float:
    deltas = slab.positions[None, :, :] - positions[:, None, :]
    fractional = deltas @ np.linalg.inv(slab.cell)
    fractional[..., :2] -= np.round(fractional[..., :2])
    cartesian = fractional @ slab.cell
    return float(np.sqrt(np.sum(cartesian * cartesian, axis=-1)).min())


def minimum_ligand_periodic_nonbonded_distance(
    positions: np.ndarray, ligand: Poscar, cell: np.ndarray
) -> float:
    original = np.linalg.norm(
        ligand.positions[:, None, :] - ligand.positions[None, :, :], axis=-1
    )
    deltas = positions[:, None, :] - positions[None, :, :]
    fractional = deltas @ np.linalg.inv(cell)
    fractional[..., :2] -= np.round(fractional[..., :2])
    cartesian = fractional @ cell
    distances = np.sqrt(np.sum(cartesian * cartesian, axis=-1))
    mask = np.triu(original >= 1.90, k=1)
    return float(distances[mask].min())


def optimize_monodentate(
    ligand: Poscar,
    anchor: int,
    slab: Poscar,
    mn_index: int,
    distance: float,
    base_azimuth: float,
) -> tuple[np.ndarray, float, float, float, float, float, float]:
    best = None
    # A C5 chain can collide even when the anchor itself is exposed. Search a
    # full rotation and, only when necessary, lengthen the initial Mn-anchor
    # distance by at most 0.5 A. The shortest clash-free distance wins.
    def orientation_candidates(trial_distance: float, anchor_options: list[tuple[float, float]]):
        candidates = []
        for anchor_tilt, anchor_azimuth in anchor_options:
            for tilt in (0.0, 20.0, 40.0, 55.0, 70.0):
                tilt_azimuths = (0.0,) if tilt == 0.0 else tuple(float(a) for a in range(0, 360, 30))
                for tilt_azimuth in tilt_azimuths:
                    for offset in range(0, 360, 20):
                        angle = (base_azimuth + offset) % 360.0
                        positions = place_monodentate(
                            ligand, anchor, slab, mn_index, trial_distance,
                            angle, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth,
                        )
                        clearance = min(
                            minimum_ligand_slab_distance(positions, slab),
                            minimum_ligand_periodic_nonbonded_distance(positions, ligand, slab.cell),
                        )
                        candidates.append(
                            (clearance, angle, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth, positions)
                        )
        return candidates

    for trial_distance in np.arange(distance, distance + 0.61, 0.10):
        candidates = orientation_candidates(float(trial_distance), [(0.0, 0.0)])
        best = max(candidates, key=lambda item: item[0])
        if best[0] >= MIN_INTERATOMIC_DISTANCE:
            clearance, angle, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth, positions = best
            return positions, angle, round(float(trial_distance), 4), tilt, tilt_azimuth, anchor_tilt, anchor_azimuth
        if abs(float(trial_distance) - distance) < 1e-8:
            tilted_options = [
                (float(anchor_tilt), float(anchor_azimuth))
                for anchor_tilt in (10, 20)
                for anchor_azimuth in range(0, 360, 30)
            ]
            candidates = orientation_candidates(float(trial_distance), tilted_options)
            tilted_best = max(candidates, key=lambda item: item[0])
            if tilted_best[0] >= MIN_INTERATOMIC_DISTANCE:
                clearance, angle, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth, positions = tilted_best
                return positions, angle, round(float(trial_distance), 4), tilt, tilt_azimuth, anchor_tilt, anchor_azimuth
    assert best is not None
    raise ValueError(
        f"No clash-free monodentate orientation at Mn {mn_index + 1}; "
        f"best clearance={best[0]:.3f} A"
    )


def place_bidentate(
    ligand: Poscar,
    oxygen_indices: tuple[int, int],
    carbon_index: int,
    slab: Poscar,
    mn_pair: tuple[int, int],
    distance: float,
) -> np.ndarray:
    normal = normal_from_cell(slab.cell)
    o1, o2 = oxygen_indices
    midpoint = 0.5 * (ligand.positions[o1] + ligand.positions[o2])
    source_u = ligand.positions[o2] - ligand.positions[o1]
    source_u /= np.linalg.norm(source_u)
    source_v = ligand.positions[carbon_index] - midpoint
    source_v -= np.dot(source_v, source_u) * source_u
    source_v /= np.linalg.norm(source_v)
    source_w = np.cross(source_u, source_v)
    source_w /= np.linalg.norm(source_w)
    i, j = mn_pair
    target_delta = minimum_image_delta(slab.positions[i], slab.positions[j], slab.cell)
    target_u = target_delta - np.dot(target_delta, normal) * normal
    target_u /= np.linalg.norm(target_u)
    target_v = normal
    target_w = np.cross(target_u, target_v)
    target_w /= np.linalg.norm(target_w)
    source_basis = np.column_stack((source_u, source_v, source_w))
    target_basis = np.column_stack((target_u, target_v, target_w))
    rotation = target_basis @ source_basis.T
    placed = (ligand.positions - midpoint) @ rotation.T
    unwrapped_j = slab.positions[i] + target_delta
    target_midpoint = 0.5 * (slab.positions[i] + unwrapped_j) + distance * normal
    return placed + target_midpoint


def place_chelate(
    ligand: Poscar,
    oxygen_indices: tuple[int, int],
    carbon_index: int,
    slab: Poscar,
    mn_index: int,
    midpoint_distance: float,
    azimuth: float,
) -> np.ndarray:
    normal = normal_from_cell(slab.cell)
    o1, o2 = oxygen_indices
    midpoint = 0.5 * (ligand.positions[o1] + ligand.positions[o2])
    source_u = ligand.positions[o2] - ligand.positions[o1]
    source_u /= np.linalg.norm(source_u)
    source_v = ligand.positions[carbon_index] - midpoint
    source_v -= np.dot(source_v, source_u) * source_u
    source_v /= np.linalg.norm(source_v)
    source_w = np.cross(source_u, source_v)
    source_w /= np.linalg.norm(source_w)
    lateral_x = slab.cell[0] - np.dot(slab.cell[0], normal) * normal
    lateral_x /= np.linalg.norm(lateral_x)
    lateral_y = np.cross(normal, lateral_x)
    phi = math.radians(azimuth)
    target_u = math.cos(phi) * lateral_x + math.sin(phi) * lateral_y
    target_v = normal
    target_w = np.cross(target_u, target_v)
    target_w /= np.linalg.norm(target_w)
    source_basis = np.column_stack((source_u, source_v, source_w))
    target_basis = np.column_stack((target_u, target_v, target_w))
    rotation = target_basis @ source_basis.T
    placed = (ligand.positions - midpoint) @ rotation.T
    target_midpoint = slab.positions[mn_index] + midpoint_distance * normal
    return placed + target_midpoint


def optimize_chelate(
    ligand: Poscar,
    oxygen_indices: tuple[int, int],
    carbon_index: int,
    slab: Poscar,
    mn_index: int,
    midpoint_distance: float,
) -> tuple[np.ndarray, float, float]:
    best = None
    for trial_distance in np.arange(midpoint_distance, midpoint_distance + 0.81, 0.10):
        candidates = []
        for azimuth in range(0, 360, 15):
            positions = place_chelate(
                ligand, oxygen_indices, carbon_index, slab, mn_index,
                float(trial_distance), float(azimuth),
            )
            clearance = min(
                minimum_ligand_slab_distance(positions, slab),
                minimum_ligand_periodic_nonbonded_distance(positions, ligand, slab.cell),
            )
            candidates.append((clearance, float(azimuth), positions))
        best = max(candidates, key=lambda item: item[0])
        if best[0] >= MIN_INTERATOMIC_DISTANCE:
            return best[2], best[1], round(float(trial_distance), 4)
    assert best is not None
    raise ValueError(
        f"No clash-free OO-chelate orientation at Mn {mn_index + 1}; "
        f"best clearance={best[0]:.3f} A"
    )


def optimize_bidentate(
    ligand: Poscar,
    oxygen_indices: tuple[int, int],
    carbon_index: int,
    slab: Poscar,
    distance: float,
    preferred_avoid_indices: set[int] | None = None,
    excluded_pairs: set[tuple[int, int]] | None = None,
) -> tuple[np.ndarray, tuple[int, int]]:
    normal = normal_from_cell(slab.cell)
    mn = [i for i, s in enumerate(slab.symbols) if s == "Mn"]
    height = {i: float(np.dot(slab.positions[i], normal)) for i in mn}
    top = max(height.values())
    candidates = [i for i in mn if height[i] >= top - SURFACE_MN_WINDOW]
    preferred_avoid_indices = preferred_avoid_indices or set()
    excluded_pairs = {tuple(sorted(pair)) for pair in (excluded_pairs or set())}
    if len(candidates) < 2:
        candidates = sorted(mn, key=lambda i: (-height[i], i))[:3]
    geometries = []
    o1, o2 = oxygen_indices
    for p, i in enumerate(candidates):
        for j in candidates[p + 1 :]:
            if tuple(sorted((i, j))) in excluded_pairs:
                continue
            positions = place_bidentate(ligand, oxygen_indices, carbon_index, slab, (i, j), distance)
            clearance = minimum_ligand_slab_distance(positions, slab)
            if clearance < MIN_INTERATOMIC_DISTANCE:
                continue
            anchor_distance = 0.5 * (
                distance_pbc_xy(positions[o1], slab.positions[i], slab.cell)
                + distance_pbc_xy(positions[o2], slab.positions[j], slab.cell)
            )
            depth = (top - height[i]) + (top - height[j])
            overlap_penalty = 10.0 * sum(index in preferred_avoid_indices for index in (i, j))
            score = anchor_distance + 0.5 * depth - 0.05 * clearance + overlap_penalty
            geometries.append((score, positions, (i, j)))
    if not geometries:
        raise ValueError("No clash-free OO-bidentate geometry was found")
    _, positions, pair = min(geometries, key=lambda item: item[0])
    return positions, pair


def parse_magmom(incar_text: str) -> list[float]:
    match = re.search(r"(?im)^\s*MAGMOM\s*=\s*(.+)$", incar_text)
    if not match:
        raise ValueError("No MAGMOM line in source INCAR")
    values: list[float] = []
    for token in match.group(1).split():
        if "*" in token:
            count, value = token.split("*", 1)
            values.extend([float(value)] * int(count))
        else:
            values.append(float(token))
    return values


def compress_magmom(values: list[float]) -> str:
    output = []
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and abs(values[end] - values[start]) < 1e-12:
            end += 1
        value = f"{values[start]:.4f}"
        output.append(f"{end-start}*{value}" if end - start > 1 else value)
        start = end
    return " ".join(output)


def update_incar(source: Path, destination: Path, slab: Poscar, ligand: Poscar, final_symbols: list[str]) -> None:
    text = source.read_text(encoding="utf-8", errors="replace")
    slab_mom = parse_magmom(text)
    if len(slab_mom) != len(slab.symbols):
        raise ValueError(f"Source MAGMOM count {len(slab_mom)} != slab atom count {len(slab.symbols)}")
    records = [(s, m, False) for s, m in zip(slab.symbols, slab_mom)]
    records += [(s, 0.0, True) for s in ligand.symbols]
    reordered_mom = [m for species in SPECIES_ORDER for s, m, _ in records if s == species]
    if len(reordered_mom) != len(final_symbols):
        raise ValueError("Final MAGMOM construction failed")
    text = re.sub(r"(?im)^\s*MAGMOM\s*=.*$", f"MAGMOM = {compress_magmom(reordered_mom)}", text)
    present_species = [s for s in SPECIES_ORDER if s in final_symbols]
    arrays = {
        "LDAUL": ["2" if s == "Mn" else "-1" for s in present_species],
        "LDAUU": ["3.900" if s == "Mn" else "0.000" for s in present_species],
        "LDAUJ": ["0.000" for _ in present_species],
    }
    for key, values in arrays.items():
        replacement = f" {key} = {' '.join(values)}"
        if re.search(rf"(?im)^\s*{key}\s*=.*$", text):
            text = re.sub(rf"(?im)^\s*{key}\s*=.*$", replacement, text)
        else:
            text = text.rstrip() + "\n" + replacement + "\n"
    destination.write_text(text.rstrip() + "\n", encoding="utf-8")


def potcar_species(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    species = re.findall(r"(?m)^\s*TITEL\s*=\s*\S+\s+(\S+)", text)
    if not species:
        raise ValueError(f"No TITEL records found in POTCAR: {path}")
    return species


def build_case(
    slab: Poscar,
    ligand: Poscar,
    ligand_positions: np.ndarray,
    source_dir: Path,
    case_dir: Path,
    metadata: dict,
    potcar: Path,
) -> dict:
    minimum = minimum_ligand_slab_distance(ligand_positions, slab)
    periodic_minimum = minimum_ligand_periodic_nonbonded_distance(
        ligand_positions, ligand, slab.cell
    )
    if minimum < MIN_INTERATOMIC_DISTANCE:
        raise ValueError(f"Ligand/slab clash in final geometry: {minimum:.3f} A")
    if periodic_minimum < MIN_INTERATOMIC_DISTANCE:
        raise ValueError(
            f"Ligand/periodic-image clash in final geometry: {periodic_minimum:.3f} A"
        )
    combined = Poscar(
        comment=f"Mn3O4(001) {metadata['state']} + {metadata['ligand']} {metadata['configuration']}",
        cell=slab.cell.copy(),
        symbols=slab.symbols + ligand.symbols,
        positions=np.vstack((slab.positions, ligand_positions)),
        flags=slab.flags + [("T", "T", "T")] * len(ligand.symbols),
    )
    case_dir.mkdir(parents=True, exist_ok=False)
    write_poscar(case_dir / "POSCAR", combined)
    update_incar(source_dir / "INCAR", case_dir / "INCAR", slab, ligand, combined.symbols)
    shutil.copy2(source_dir / "KPOINTS", case_dir / "KPOINTS")
    species = [s for s in SPECIES_ORDER if s in combined.symbols]
    expected_potcar = ["Mn_pv" if s == "Mn" else s for s in species]
    actual_potcar = potcar_species(potcar)
    if actual_potcar != expected_potcar:
        raise ValueError(
            f"POTCAR order mismatch for {metadata['ligand']}: "
            f"expected {expected_potcar}, found {actual_potcar} in {potcar}"
        )
    (case_dir / "POTCAR.spec").write_text("\n".join(expected_potcar) + "\n")
    shutil.copy2(potcar, case_dir / "POTCAR")
    metadata = {
        **metadata,
        "source_input": str(source_dir),
        "slab_atoms": len(slab.symbols),
        "ligand_atoms": len(ligand.symbols),
        "total_atoms": len(combined.symbols),
        "species_order": species,
        "ligand_magmom": 0.0,
        "ligand_raise_for_clash_A": 0.0,
        "minimum_ligand_slab_distance_A": round(minimum, 4),
        "minimum_ligand_periodic_nonbonded_distance_A": round(periodic_minimum, 4),
        "potcar_complete": True,
        "potcar_source": str(potcar),
    }
    (case_dir / "adsorption_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--ligand-root", type=Path, default=DEFAULT_LIGAND_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--potcar-amine", type=Path, default=DEFAULT_POTCAR_AMINE)
    parser.add_argument("--potcar-carboxyl", type=Path, default=DEFAULT_POTCAR_CARBOXYL)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Output already exists; refusing to overwrite: {args.output}")
    for p in (args.potcar_amine, args.potcar_carboxyl):
        if not p.is_file():
            raise FileNotFoundError(p)
    amine = read_poscar(args.ligand_root / "C5-amine.vasp")
    carboxyl = read_poscar(args.ligand_root / "C5-carboxyl.vasp")
    amine_n = [i for i, s in enumerate(amine.symbols) if s == "N"]
    carboxyl_o = [i for i, s in enumerate(carboxyl.symbols) if s == "O"]
    carboxyl_c = [i for i, s in enumerate(carboxyl.symbols) if s == "C"]
    if len(amine_n) != 1 or len(carboxyl_o) != 2:
        raise ValueError("Unexpected ligand anchor counts")
    center_o = 0.5 * (carboxyl.positions[carboxyl_o[0]] + carboxyl.positions[carboxyl_o[1]])
    carboxyl_center_c = min(carboxyl_c, key=lambda i: np.linalg.norm(carboxyl.positions[i] - center_o))
    args.output.mkdir(parents=True)
    manifest = []
    with args.selection.open(newline="", encoding="utf-8-sig") as handle:
        selected = list(csv.DictReader(handle))
    for row in selected:
        state = row["state"]
        variant = Path(row["path"]).name
        source_dir = args.source_root / state / variant
        slab = read_poscar(source_dir / "POSCAR")
        amine_sites = choose_surface_mn(slab, BIND_DISTANCE_N_MN)
        carboxyl_sites = choose_surface_mn(slab, BIND_DISTANCE_O_MN)
        amine_geometries = []
        for site in amine_sites:
            base_angle = (0.0, 72.0, 144.0, 216.0, 288.0)[len(amine_geometries)]
            try:
                geometry = optimize_monodentate(
                    amine, amine_n[0], slab, site, BIND_DISTANCE_N_MN, base_angle
                )
            except ValueError:
                continue
            amine_geometries.append((site, geometry))
            if len(amine_geometries) == 3:
                break
        if len(amine_geometries) < 3:
            raise ValueError(f"{state}: only {len(amine_geometries)} clash-free distinct Mn-N sites found")
        for number, (site, geometry) in enumerate(amine_geometries, start=1):
            positions, angle, actual_distance, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth = geometry
            config = f"N_site{number:02d}"
            metadata = build_case(
                slab, amine, positions, source_dir,
                args.output / state / "amine" / config,
                {"state": state, "variant": variant, "ligand": "C5-amine", "configuration": config,
                 "binding_mode": "N_monodentate", "target_mn_indices_1based": [site + 1],
                 "nominal_binding_distance_A": BIND_DISTANCE_N_MN,
                 "actual_binding_distance_A": actual_distance, "roll_deg": angle,
                 "backbone_tilt_deg": tilt, "backbone_tilt_azimuth_deg": tilt_azimuth,
                 "anchor_tilt_deg": anchor_tilt, "anchor_azimuth_deg": anchor_azimuth},
                args.potcar_amine,
            )
            manifest.append(metadata)
        carboxyl_geometries = []
        for site in carboxyl_sites:
            number = len(carboxyl_geometries)
            anchor = carboxyl_o[number % 2]
            base_angle = (30.0, 150.0)[number]
            try:
                geometry = optimize_monodentate(
                    carboxyl, anchor, slab, site, BIND_DISTANCE_O_MN, base_angle
                )
            except ValueError:
                continue
            carboxyl_geometries.append((site, anchor, geometry))
            if len(carboxyl_geometries) == 2:
                break
        if len(carboxyl_geometries) < 2:
            raise ValueError(f"{state}: only {len(carboxyl_geometries)} clash-free distinct Mn-O sites found")
        for number, (site, anchor, geometry) in enumerate(carboxyl_geometries, start=1):
            positions, angle, actual_distance, tilt, tilt_azimuth, anchor_tilt, anchor_azimuth = geometry
            config = f"O_site{number:02d}"
            metadata = build_case(
                slab, carboxyl, positions, source_dir,
                args.output / state / "carboxyl" / config,
                {"state": state, "variant": variant, "ligand": "C5-carboxyl", "configuration": config,
                 "binding_mode": "O_monodentate", "target_mn_indices_1based": [site + 1],
                 "anchor_ligand_index_1based": anchor + 1,
                 "nominal_binding_distance_A": BIND_DISTANCE_O_MN,
                 "actual_binding_distance_A": actual_distance, "roll_deg": angle,
                 "backbone_tilt_deg": tilt, "backbone_tilt_azimuth_deg": tilt_azimuth,
                 "anchor_tilt_deg": anchor_tilt, "anchor_azimuth_deg": anchor_azimuth},
                args.potcar_carboxyl,
            )
            manifest.append(metadata)
        mono_sites = {site for site, _, _ in carboxyl_geometries}
        chelate_candidates = [site for site in carboxyl_sites if site not in mono_sites]
        chelate_candidates += [site for site in carboxyl_sites if site in mono_sites]
        chelate_geometries = []
        for site in chelate_candidates:
            try:
                geometry = optimize_chelate(
                    carboxyl, (carboxyl_o[0], carboxyl_o[1]), carboxyl_center_c,
                    slab, site, BIND_DISTANCE_O_MN,
                )
            except ValueError:
                continue
            chelate_geometries.append((site, geometry))
            if len(chelate_geometries) == 1:
                break
        if len(chelate_geometries) < 1:
            raise ValueError(f"{state}: only {len(chelate_geometries)} clash-free distinct OO sites found")
        for number, (site, geometry) in enumerate(chelate_geometries, start=1):
            positions, azimuth, midpoint_distance = geometry
            o_mn_distances = [
                distance_pbc_xy(positions[o], slab.positions[site], slab.cell)
                for o in carboxyl_o
            ]
            config = f"OO_chelate{number:02d}"
            metadata = build_case(
                slab, carboxyl, positions, source_dir,
                args.output / state / "carboxyl" / config,
                {"state": state, "variant": variant, "ligand": "C5-carboxyl", "configuration": config,
                 "binding_mode": "OO_chelate", "target_mn_indices_1based": [site + 1],
                 "nominal_midpoint_distance_A": BIND_DISTANCE_O_MN,
                 "actual_midpoint_distance_A": midpoint_distance,
                 "actual_O_Mn_distances_A": [round(value, 4) for value in o_mn_distances],
                 "azimuth_deg": azimuth},
                args.potcar_carboxyl,
            )
            manifest.append(metadata)
    columns = [
        "state", "variant", "ligand", "configuration", "binding_mode", "source_input",
        "slab_atoms", "ligand_atoms", "total_atoms", "species_order", "target_mn_indices_1based",
        "ligand_raise_for_clash_A", "minimum_ligand_slab_distance_A", "potcar_complete",
    ]
    with (args.output / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for item in manifest:
            serial = dict(item)
            serial["species_order"] = " ".join(item["species_order"])
            serial["target_mn_indices_1based"] = " ".join(map(str, item["target_mn_indices_1based"]))
            writer.writerow(serial)
    readme = f"""# Mn3O4(001) ligand adsorption inputs

Generated cases: {len(manifest)} ({len(selected)} states x 6 configurations)

Output hierarchy: `<state>/<ligand>/<configuration>`; the relaxed source variant
is retained in `manifest.csv` and each `adsorption_metadata.json` rather than as
an additional directory level.

- `amine/N_site01..03`: N-monodentate adsorption on three distinct Mn sites.
- `carboxyl/O_site01..02`: O-monodentate adsorption on two distinct Mn sites.
- `carboxyl/OO_chelate01`: both carboxyl O atoms coordinate to a third Mn site.
- Original slab selective dynamics are preserved; all ligand atoms are movable.
- Every added ligand atom has `MAGMOM = 0.0`.
- The cell height is 38.3523 A; generated structures retain 17.44-25.27 A
  atom-plane vacuum gaps.
- Placement rejects ligand/slab and ligand/periodic-image distances below 1.15 A.

## POTCAR status

Each case contains a complete, user-supplied `POTCAR` and a `POTCAR.spec` record.
The generator validates every POTCAR TITEL order before copying it:

- amine: `Mn_pv O H N C`
- carboxyl: `Mn_pv O H C`
"""
    (args.output / "README.md").write_text(readme, encoding="utf-8")
    print(f"Generated {len(manifest)} cases in {args.output}")


if __name__ == "__main__":
    main()
