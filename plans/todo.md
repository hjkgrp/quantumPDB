# QuantumPDB implementation to-dos

Tracks work planned/in progress on this branch. Each chapter is one
self-contained change (its own commit).

## 1. `force_include_residues`: force-include specific protein residues

**Plan:** In the spirit of `additional_ligands`, let users force specific
protein residues (identified as `RESNAME_CHAINID`, e.g. `HIS_A123`) into the
QM cluster even if sphere growth (distance cutoff + Voronoi adjacency) never
reaches them — they may lie beyond the first coordination sphere. Added
residues must be properly extracted and capped like any other cluster
residue, and protected from `max_atom_count` pruning. Configurable both
globally in YAML (`force_include_residues: [...]`, applied to every PDB) and
per-PDB via a new CSV column (takes priority over YAML when present),
mirroring the existing `center`/`center_residues` precedence.

Full design: `plans/force_include_residues.md` (copied from the approved plan
file; see its "Implementation Summary" section for what actually shipped,
including two deviations from the original draft below).

**Status: implemented, tests passing, not yet committed.**

- [x] `qp/cluster/spheres.py`: add `add_force_include_residues()`, extend
      `prune_atoms()` with `protected_residues`, wire both into
      `extract_clusters()` (new `force_include_residues=[]` param).
- [x] `qp/structure/setup.py`: add a standalone `get_force_include_residues()`
      CSV-column reader (kept `get_centers()` untouched, per review — see
      Implementation Summary); extend `parse_input()` to resolve per-PDB
      `force_include_residues` from CSV or YAML.
- [x] `qp/cli.py`: read `force_include_residues` from config, thread it through
      `parse_input()` (in both `run()` and `analyze()`), pop it per-PDB
      alongside `center_residues` (incl. in the error/skip branches), pass
      into `extract_clusters()`.
- [x] `config.yaml` / `qp/tests/samples.yaml`: add `force_include_residues: []`
      template entry next to `additional_ligands`.
- [x] `docs/configuration.rst`: document the new parameter.
- [x] `docs/input_formats.rst`: document the new optional CSV column and its
      dash-separated multi-residue syntax, plus the YAML list syntax.
- [x] `qp/tests/test_cluster.py`: add `test_force_include_residues` (inclusion
      beyond default spheres, capping, survives `max_atom_count` pruning),
      `test_force_include_residues_multiple` (two residues, two chains, one
      call), and `test_force_include_residues_already_present` (force-including
      a residue already in the cluster doesn't write it twice).
- [x] Renamed `additional_residues` → `force_include_residues` everywhere
      (config key, CSV column, all internal identifiers, docs headings/
      anchors, this plan file's own name) — repo-wide grep confirms no
      leftover references.
- [x] Verify: ran all three new tests from the `my_enerzyme_dev` conda env
      (has Bio/sklearn/scipy) — all pass. Ran the full `test_cluster.py`
      suite too: 11 pre-existing failures, confirmed (via `git stash`/`pop`
      against the original commit) to predate this change — they pass a raw
      `["FE", "FE2"]` list as `center_residue` where the current API needs a
      `CenterResidue` object. Unrelated to this feature; not fixed here.

## 2. `additional_ligands` CSV column (deferred, separate commit)

**Plan:** Extend `additional_ligands` — currently YAML-only, computed once
in `qp/cli.py` and applied identically to every PDB in a batch — with the
same per-PDB CSV override mechanism built for `force_include_residues` above
(a standalone `get_additional_ligands()` reader, following the
`get_force_include_residues()` precedent, plus the matching `parse_input()`
resolution logic). Lets a batch mix structures that each need different
extra ligand/cofactor resnames kept, instead of one global list for the
whole run.

- [ ] `qp/structure/setup.py`: add a standalone `get_additional_ligands()`
      CSV-column reader and resolve it in `parse_input()`, same
      CSV-takes-priority-over-YAML precedence as `force_include_residues`.
- [ ] `qp/cli.py`: move `ligands` resolution from "compute once before the
      loop" to "pop per-PDB" (like `center_residues`/`force_include_residues`),
      including the error/skip branches.
- [ ] `config.yaml` / docs: document the new CSV column.
- [ ] Tests: cover per-PDB ligand override via CSV.
