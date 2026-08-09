# Input File Agent Memory

Last updated: 2026-08-05

## Current Repository Facts

- Workspace root: `/home/sowon-desktop/EnergeticSpanModel`
- This directory is not currently detected as a git repository by `git status`.
- `.agents/` and `.codex/` existed before this setup, but they had no instruction files inside.
- Root `AGENTS.md` now points future agents to this folder for VASP input-generation work.

## Reusable Code

- Shared generator engine: `011_facet/generate_mn3o4_state_inputs.py`
- 001 facet wrapper: `001_facet/generate_mn3o4_state_inputs.py`
- 011 34-layer wrapper: `011_facet_34/generate_mn3o4_state_inputs.py`
- Older direct MAGMOM updater example: `011_facet/12layer_3rd_magmom/add_magmom3.py`
- Other older helpers exist under `001_facet/vac30_layer12_magmom/` and `011_facet/12layer_magmom/`.

## Generated Input Directories Observed

- `001_facet/generated_state_inputs_001/`
- `011_facet/generated_state_inputs_2/`
- `011_facet_34/generated_state_inputs_011_34layer/`

Generated variant folders usually contain:

- `POSCAR`
- `INCAR`
- `KPOINTS`
- `POTCAR`
- `state_metadata.json`

Some generated roots also include `manifest.csv` and `README.md`.

## MAGMOM Conventions Observed

001 facet base MAGMOM from `001_facet/magmom_test/initial/INCAR`:

```text
4*3.8000 4*-3.8000 4*4.4
```

In wrapper form:

```python
BASE_MN_MAGMOM = {
    0: 3.8, 1: 3.8, 2: 3.8, 3: 3.8,
    4: -3.8, 5: -3.8, 6: -3.8, 7: -3.8,
    8: 4.4, 9: 4.4, 10: 4.4, 11: 4.4,
}
```

011 12-layer base MAGMOM in shared engine:

```python
BASE_MN_MAGMOM = {
    0: 3.8, 1: 3.8, 2: -3.8, 3: -3.8,
    4: 3.8, 5: 3.8, 6: -3.8, 7: -3.8,
    8: 4.4, 9: 4.4, 10: 4.4, 11: 4.4,
}
```

Older direct updater example for third-layer 011 cases:

```text
2*3.8000 2*-3.8000 2*3.8000 2*-3.8000 3*4.4 0.0 4.4
```

## Durable Assumptions

- Mn slab atoms keep nonzero MAGMOM.
- Added/generated non-Mn atoms usually get `0.0000`.
- Generated INCAR should remove previous `MAGMOM` and `NCORE` lines before appending fresh values.
- `NCORE = 8` has been used in older direct updaters.
- The generator compresses repeated MAGMOM values, for example `4*3.8000`.

## How To Update This Memory

Add a new dated entry under "Change Log" whenever:

- A new generator or wrapper is created.
- A base MAGMOM rule changes.
- A new facet, layer count, or reference directory is introduced.
- Output directory naming changes.
- A validation command is run and its result matters.

## Change Log

### 2026-08-05

- Created `.agents/input_file_agent/` as the persistent folder for input-file generation workflow notes.
- Created root `AGENTS.md` so future agents know where to read VASP input-generation instructions.
- Added second-pass 001 CONTCAR-to-input workflow:
  - Source CONTCAR folder: `generated_001_contcar/`
  - Reference bottom-layer structure: `generated_001_contcar/contcar_s000_01.vasp`
  - Magnetization source: `generated_001_contcar/OUTCAR_322`
  - Generator script: `generated_001_contcar/make_initial_inputs.py`
  - Output inputs: `generated_001_contcar_initial_inputs/`
  - Rule: replace/fix original atoms below `z < 14.0 Angstrom` with positions from `contcar_s000_01.vasp`.
  - Rule: assign Mn MAGMOM by position-matching each state Mn to Mn sites in `contcar_s322_02.vasp` with `OUTCAR_322` magnetization, then map to `+/-3.8`, `+/-4.4`, or `0.0`.
  - Generated 26 one-variant input sets preserving filename variants, e.g. `contcar_s102_02.vasp` -> `S102/variant_02`.
  - Validation: MAGMOM count matched atom count for all 26 sets; bottom fixed positions matched the reference within numerical precision.
