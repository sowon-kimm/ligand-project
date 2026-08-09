from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


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


def main() -> None:
    generator = load_generator()

    generator.HERE = HERE
    generator.START_DIR = HERE / "magmom_test" / "initial"
    generator.FINAL_DIR = HERE / "vac30_layer12_magmom" / "olaoxioxooxioxo"
    generator.OUTPUT_DIR = HERE / "generated_state_inputs"

    # From 001_facet/magmom_test/initial/INCAR:
    # MAGMOM = 4*3.8000 4*-3.8000 4*4.4 19*0.0000
    generator.BASE_MN_MAGMOM = {
        0: 3.8,
        1: 3.8,
        2: 3.8,
        3: 3.8,
        4: -3.8,
        5: -3.8,
        6: -3.8,
        7: -3.8,
        8: 4.4,
        9: 4.4,
        10: 4.4,
        11: 4.4,
    }

    generator.main()


if __name__ == "__main__":
    main()
