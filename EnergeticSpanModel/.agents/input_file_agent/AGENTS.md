# Input File Agent

Purpose: maintain and extend the Mn3O4 VASP input-file generation workflow in this repository.

When the user asks to make input files, set INCAR MAGMOM, generate Mn3O4 state variants, or update VASP templates:

1. Read this folder first.
2. Prefer the existing generator engine in `011_facet/generate_mn3o4_state_inputs.py`.
3. Use facet-specific wrapper scripts instead of duplicating the full engine.
4. Preserve existing generated data unless the user explicitly asks to regenerate or overwrite.
5. Record durable decisions and new paths in `.agents/input_file_agent/MEMORY.md`.

Important local pattern:

- `011_facet/generate_mn3o4_state_inputs.py` is the shared engine.
- `001_facet/generate_mn3o4_state_inputs.py` imports the shared engine and overrides paths/MAGMOM.
- `011_facet_34/generate_mn3o4_state_inputs.py` imports the shared engine and derives base Mn MAGMOM from `OUTCAR`.

Default safety:

- Before overwriting generated input folders, inspect whether they already contain user data.
- Do not delete VASP output files such as `OUTCAR`, `vasprun.xml`, `WAVECAR`, or `CHGCAR` unless explicitly requested.
- Check `INCAR` and `POSCAR` atom counts when changing MAGMOM.
