# Input File Generation Runbook

Use this when asked to generate or update VASP input files.

## 1. Identify The Case

Find the facet/layer target and confirm these paths:

- start/reference structure directory
- final/reference structure directory
- output directory
- whether the final structure file is `POSCAR` or `CONTCAR`
- source directories for template files `INCAR`, `KPOINTS`, `POTCAR`

Prefer wrappers that override the shared engine instead of copying the engine.

## 2. Existing Generators

Run from repository root unless a script expects another working directory.

```bash
python 001_facet/generate_mn3o4_state_inputs.py
python 011_facet/generate_mn3o4_state_inputs.py
python 011_facet_34/generate_mn3o4_state_inputs.py
```

Before running, open the script and verify `START_DIR`, `FINAL_DIR`, `OUTPUT_DIR`, and MAGMOM rules. Some historical generated directory names may differ from current script defaults.

## 3. MAGMOM Checks

After generation, inspect several representative INCAR files:

```bash
rg -n "^MAGMOM|^NCORE" 001_facet/generated_state_inputs_001 011_facet/generated_state_inputs_2 011_facet_34/generated_state_inputs_011_34layer
```

For a specific variant:

```bash
sed -n '1,80p' path/to/variant/INCAR
```

If modifying MAGMOM manually:

1. Read atom counts from POSCAR.
2. Remove existing `MAGMOM` and `NCORE` lines from INCAR.
3. Append the new `MAGMOM = ...` line.
4. Append `NCORE = 8` only if that is still desired for the target workflow.

## 4. Validation

Useful checks:

```bash
python 001_facet/generate_mn3o4_state_inputs.py
python 011_facet/generate_mn3o4_state_inputs.py
python 011_facet_34/generate_mn3o4_state_inputs.py
rg -n "^MAGMOM" 001_facet/generated_state_inputs_001 011_facet/generated_state_inputs_2 011_facet_34/generated_state_inputs_011_34layer
find 001_facet/generated_state_inputs_001 011_facet/generated_state_inputs_2 011_facet_34/generated_state_inputs_011_34layer -name state_metadata.json | wc -l
```

If dependencies are missing, likely packages are:

- `ase`
- `numpy`
- `scipy`

## 5. When Adding A New Facet Wrapper

Create a small wrapper script beside the new facet directory. It should:

1. Import `011_facet/generate_mn3o4_state_inputs.py` with `importlib.util`.
2. Override `HERE`, `START_DIR`, `FINAL_DIR`, `OUTPUT_DIR`.
3. Override `START_STRUCTURE_FILE` or `FINAL_STRUCTURE_FILE` if needed.
4. Set `TEMPLATE_SOURCE_DIRS` if templates should come from multiple places.
5. Set `BASE_MN_MAGMOM` explicitly or derive it from `OUTCAR`.
6. Call `generator.main()`.

After adding the wrapper, update `.agents/input_file_agent/MEMORY.md`.
