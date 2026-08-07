# Add `force_include_residues`: force-include specific protein residues in the QM cluster

## Context

`qp/cluster/spheres.py` already lets users force specific *ligand/cofactor*
residues into a cluster via `additional_ligands` (config) / `ligands` (code
param): matching resnames are protected from `prune_atoms` and, when
`include_ligands == 0`, are allowed past the first sphere. There is no
equivalent for ordinary **protein residues** — a user cannot say "always
include chain A residue 245, even though it's outside the coordination
spheres that were actually grown, and even if `max_atom_count` would
otherwise prune it." This matters because the sphere-growth logic
(`get_next_neighbors`) only reaches residues via a distance cutoff (sphere 1)
or Voronoi-graph adjacency to the previous sphere (sphere 2+) — a residue
that is spatially/functionally relevant but not graph-connected in that way
(e.g. a distal catalytic or second-shell residue) can never be picked up, no
matter how `number_of_spheres` is tuned, without dragging in everything
between it and the center.

This plan adds `force_include_residues`, in the spirit of `additional_ligands`,
to force-include named protein residues regardless of how they were (or
weren't) reached by sphere growth, ensure they are properly capped, and
protect them from `max_atom_count` pruning. (A second feature — restricting
the cluster to only center-interacting residues — was discussed and
explicitly dropped from this plan; not implemented here.)

Decisions already confirmed with the user:
- Residues are specified in the same `RESNAME_CHAINID` format `spheres.py`
  already uses internally (`make_res_key`, e.g. `HIS_A123`), the same style
  as `CenterResidue`'s strict-mode tokens (e.g. `CU_A357`).
- Supported both as a global YAML list (`force_include_residues: [...]`,
  applied identically to every PDB in a run — like `additional_ligands`
  today) **and** as a per-PDB CSV column (like the existing `center`
  column), since chain+residue-number identifiers are structure-specific.
  When the CSV column is present, it takes priority over the YAML list, for
  every PDB in that run — mirroring the documented `center`/`center_residues`
  precedence exactly.
- Explicitly-requested residues are always protected from `max_atom_count`
  pruning (same guarantee `additional_ligands` already gives).

## Design

### `qp/cluster/spheres.py`

Add one new function, `add_force_include_residues`, placed right after
`complete_oligomer` (same category of step: supplement `residues`/`spheres`
before pruning/capping/output). Modeled directly on `complete_oligomer`'s
shape — iterate `model`'s chains/residues, match by `make_res_key`, mutate
`residues`/`spheres` in place, print a line per outcome — reusing
`make_res_key` rather than introducing new parsing:

```python
def add_force_include_residues(model, residues, spheres, force_include_residues):
    """Force-include specific protein residues, even beyond the grown spheres.

    Matches each entry (``'RESNAME_CHAINID'``, e.g. ``'HIS_A123'``) against
    every residue in the model. Matches not already part of the cluster are
    added to ``residues`` and to the outermost sphere. All matches (whether
    newly added or already present) are returned so callers can protect them
    from later pruning.
    ...
    Returns
    -------
    set
        Residues matching ``force_include_residues`` (new or pre-existing).
    """
    requested = set(force_include_residues)
    matched = set()
    if not requested:
        return matched
    for chain in model:
        for res in chain.get_unpacked_list():
            res_key = make_res_key(res)
            if res_key in requested:
                matched.add(res)
                requested.discard(res_key)
                if res not in residues:
                    residues.add(res)
                    spheres[-1].add(res)
                    print(f"> {res_key} added to cluster via force_include_residues")
    for res_key in requested:
        print(f"> WARNING: force_include_residues entry {res_key!r} was not found in the structure")
    return matched
```

Extend `prune_atoms` with one new parameter so it can protect these
residues by identity (resname alone isn't enough — unlike ligands, an added
protein residue commonly shares its resname, e.g. `HIS`, with many
unrelated residues already in the cluster):

```python
def prune_atoms(center, residues, spheres, max_atom_count, ligands, protected_residues=frozenset()):
    ...
    for res in sorted(residues, key=dist, reverse=True):
        if res.get_resname() not in ligands and res not in protected_residues:
            prune.add(res)
            ...
```

In `extract_clusters`, add an `force_include_residues=[]` parameter (same
default style as `ligands`), documented in the docstring next to the
existing `ligands` entry. Inside the per-center loop, right after the
existing `complete_oligomer(...)` call:

```python
added_residues = add_force_include_residues(model, residues, spheres, force_include_residues)
...
if max_atom_count is not None:
    prune_atoms(c, residues, spheres, max_atom_count, ligands, added_residues)
```

No changes needed to `cap_chains`, `compute_charge`, `count_residues`,
`write_pdbs`, or `struct_to_file`: they all operate on `residues`/`spheres`
by membership, not by how a residue got there, so newly-added residues are
picked up, capped (real chain-break detection against the original,
uncut chain), charge-counted, and written to the outermost sphere's PDB /
combined `.xyz`/`.pdb` automatically. Also update the module's top-of-file
usage docstring example to mention the new parameter, mirroring the
existing `ligands=["AKG"]` line.

### `qp/structure/setup.py`

Reuse `get_centers` for the new CSV column instead of writing a near-duplicate
reader: give it a `column` parameter (default `'center'`, so the existing
call site and its behavior are unchanged) and a `keep_blank` flag (CSV rows
with no `center` value are currently skipped outright — fine for a required
column — but `force_include_residues` is optional per-row, so blank cells must
still produce a placeholder to keep row alignment with `pdb_all`):

```python
def get_centers(input_path, column='center', keep_blank=False):
    ...
    for row in reader:
        value = row.get(column, None)
        if value:
            centers.append(value)
        elif keep_blank:
            centers.append("")
```

Extend `parse_input` to accept `force_include_yaml_residues` and resolve the
final per-PDB list the same way it already resolves `center_residues`,
so `cli.py` doesn't need to know about CSV-vs-YAML precedence at all:

```python
def parse_input(input, output, center_yaml_residues, force_include_yaml_residues=None):
    ...
    force_include_csv_residues = get_centers(input, column='force_include_residues', keep_blank=True)
    ...
    if force_include_csv_residues:
        print("> Using force-included residues from the input csv\n")
        force_include_residues = [[t for t in row.split("-") if t] for row in force_include_csv_residues]
    else:
        force_include_residues = [list(force_include_yaml_residues or []) for _ in pdb_all]

    return pdb_all, center_residues, force_include_residues
```

(Dash as the multi-residue delimiter mirrors the existing merged-center CSV
syntax, e.g. `CU_A357-CU_A358`.)

### `qp/cli.py`

- `center_yaml_residues = config_data.get('center_residues', [])` gains a
  sibling `force_include_yaml_residues = config_data.get('force_include_residues', [])`,
  and the `parse_input` call becomes
  `pdb_all, center_residues, force_include_residues_all = setup.parse_input(input, output, center_yaml_residues, force_include_yaml_residues)`.
- Inside the per-PDB loop, mirror the existing `center_residues.pop(0)`
  pattern exactly: pop `force_include_residues_all` alongside `center_residues`
  in **every** branch that currently does `if center_residues:
  center_residues.pop(0)` before `continue` (the oversized/`ValueError`/`IOError`
  skip branches), so the two lists stay index-aligned as PDBs are skipped.
- Where `center_residue = CenterResidue(center_residues.pop(0), ...)` is
  built, add `force_include_residues = force_include_residues_all.pop(0) if force_include_residues_all else []`.
- Pass `force_include_residues=force_include_residues` into the
  `spheres.extract_clusters(...)` call (as a keyword arg, next to the
  existing `cluster_name_template=...` keyword).

`additional_ligands` (`ligands`) is untouched — it stays YAML-only, computed
once before the loop, exactly as today.

### Docs & config template

- `config.yaml`: add `force_include_residues: []` under the cluster model
  parameters, next to `additional_ligands`, with a comment describing the
  `RESNAME_CHAINID` format.
- `docs/configuration.rst`: add a table row for `force_include_residues` next
  to `additional_ligands`.
- `docs/input_formats.rst`: document the new optional CSV column next to
  `center`/`oxidation`/`multiplicity`, including the dash-separated
  multi-residue syntax.
- `qp/tests/samples.yaml`: add `force_include_residues: []` next to
  `additional_ligands: []` for template consistency.

### Tests

Add `test_force_include_residues` to `qp/tests/test_cluster.py`, following the
shape of `test_prune_atoms`/`test_cap_heavy`: pick a residue from one of the
existing sample structures that is *not* included by default, run
`extract_clusters(..., force_include_residues=["<RESNAME>_<CHAINID>"], capping=1)`,
and assert the residue's PDB record appears in the outermost sphere's
output and that it received a cap. Also verify it survives when
`max_atom_count` is set low enough that it would otherwise be pruned.

## Verification

- `python -c "import Bio"` failed in this session's active interpreter (no
  `Bio`/`sklearn`/`scipy`) — the project's actual dependencies live in a
  dedicated conda/venv env per `environment.yml`/`pyproject.toml`, not the
  base one. Run tests from that project env:
  `pytest qp/tests/test_cluster.py -k force_include_residues -q` (and the full
  `qp/tests/test_cluster.py` to check for regressions — note prior
  exploration found some tests there already call `extract_clusters` with a
  stale `center_residue` signature, so pre-existing failures unrelated to
  this change may surface; report those separately rather than folding them
  into this change).
- Manual smoke test: run `qp run -c config.yaml` (or a small test config)
  against one of the sample PDBs in `qp/tests/samples/` with
  `force_include_residues: [SOME_RESIDUE_A123]` set to a residue known to sit
  outside the default spheres, and inspect the output cluster directory to
  confirm the residue's atoms are present in the outermost sphere PDB and
  in `count.csv`, and that a cap was added where its neighbor was excluded.

## Implementation Summary

What actually shipped, and where it diverged from the draft above.

**Files touched:** `qp/cluster/spheres.py`, `qp/structure/setup.py`,
`qp/cli.py`, `config.yaml`, `qp/tests/samples.yaml`,
`docs/configuration.rst`, `docs/input_formats.rst`,
`qp/tests/test_cluster.py`. Nothing committed yet.

**Two deviations from the draft, both from review feedback during
implementation:**

1. **CSV reader stayed a standalone function.** The draft proposed
   generalizing `get_centers()` with `column`/`keep_blank` params and
   reusing it for the new CSV column. Rejected on review (unnecessary
   churn to a working, unrelated function) — implemented instead as a
   separate `get_force_include_residues()` in `qp/structure/setup.py`,
   parallel to but independent from `get_centers()`. `get_centers()` is
   untouched.
2. **YAML stayed a list, not a dash-joined string.** Briefly explored
   making the YAML value a single dash-separated string to exactly mirror
   the CSV cell format (`force_include_residues: HIS_A123-GLU_A45`).
   Reverted on review: YAML keeps native list syntax
   (`force_include_residues: [HIS_A123, GLU_A45]`), consistent with
   `center_residues`. CSV still uses a dash-separated string per cell,
   since a CSV cell can't hold a YAML list — the two formats share the
   same per-token `RESNAME_CHAINID` syntax, not the same container syntax.

**Renamed after initial implementation:** `additional_residues` →
`force_include_residues`, across the config key, CSV column, and every
internal identifier (`add_force_include_residues`,
`get_force_include_residues`, `force_include_yaml_residues`,
`force_include_csv_residues`, `force_include_residues_all`, docstrings,
docs headings/anchors, this plan file's own name). Verified with a
repo-wide grep that no `additional_residues` references remain anywhere.

**Test coverage** (`qp/tests/test_cluster.py`, all against sample `3a8g`,
center `FE`):

- `test_force_include_residues` — a residue outside the default spheres
  (`GLU_A9`) is absent by default, present and capped when force-included,
  and survives `max_atom_count` pruning aggressive enough to normally
  prune it first.
- `test_force_include_residues_multiple` — two residues on different
  chains (`GLU_A9`, `GLY_B3`) both get force-included and capped in one
  call.
- `test_force_include_residues_already_present` — force-including a
  residue that's *already* part of the default cluster (`ALA_A113`, fully
  flanked by in-cluster neighbors so it isn't capped either way) does not
  duplicate it: the same atom-line count in the combined output PDB with
  or without `force_include_residues` set. Added on request, to explicitly
  regression-test the `if res not in residues:` guard in
  `add_force_include_residues`.

**Verification run:** all three tests pass from the `my_enerzyme_dev`
conda env (`Bio`/`sklearn`/`scipy` available there; the base env does
not have them). Ran the full `test_cluster.py` suite too: 11 pre-existing
failures, confirmed via `git stash`/`git stash pop` against the
unmodified commit to predate this change entirely — those tests pass a
raw `["FE", "FE2"]` list as `center_residue` where the current API
requires a `CenterResidue` object. Unrelated to this feature; not fixed
as part of it.
