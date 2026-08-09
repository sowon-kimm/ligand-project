from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

from ase.io import read


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ENGINE_PATH = ROOT / "011_facet" / "generate_mn3o4_state_inputs.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("mn3o4_state_input_generator", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load generator from {ENGINE_PATH}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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


def base_mn_magmom_from_outcar(poscar_path: Path, outcar_path: Path) -> dict[int, float]:
    atoms = read(poscar_path)
    tots = parse_last_magnetization_tot(outcar_path)
    if len(tots) < len(atoms):
        raise RuntimeError(f"OUTCAR has {len(tots)} magnetization rows, but POSCAR has {len(atoms)} atoms")

    return {
        idx: magmom_from_tot(tots[idx])
        for idx, atom in enumerate(atoms)
        if atom.symbol == "Mn"
    }


def main() -> None:
    generator = load_generator()

    generator.HERE = HERE
    generator.START_DIR = HERE / "12layer_final"
    generator.FINAL_DIR = HERE / "34layer"
    generator.OUTPUT_DIR = HERE / "generated_state_inputs"
    generator.START_STRUCTURE_FILE = "POSCAR"
    generator.FINAL_STRUCTURE_FILE = "CONTCAR"
    generator.TEMPLATE_SOURCE_DIRS = (generator.START_DIR, generator.FINAL_DIR)
    generator.BASE_MN_MAGMOM = base_mn_magmom_from_outcar(
        generator.START_DIR / generator.START_STRUCTURE_FILE,
        generator.START_DIR / "OUTCAR",
    )

    generator.main()


if __name__ == "__main__":
    main()
