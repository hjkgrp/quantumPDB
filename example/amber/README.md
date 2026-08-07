# AMBER cluster extraction example

Extract QM coordination-sphere clusters from a user-supplied, AMBER ff14SB
prepared structure (`pdb/small_amber.pdb`) using two equivalent entry points.

## Active development
- Handling of user-supplied charge
- Point charge generation/handling
- Warnings of other nonstandard residues that users might need to check

## What is special about AMBER input

- AMBER protonation-state residue names (`HIE`/`HIP`/`HID`, `CYX`, `ASH`/`GLH`,
  `LYN`, ...) are **automatically renamed** to canonical PDB names for any local
  `.pdb` input, so BioPython's `is_aa` and downstream steps work. The untouched
  original is preserved as `output/small_amber/small_amber_original.pdb`.
- The structure is already modeled and protonated, so **Modeller and Protoss are
  skipped**. Capping, charge counting, and residue counting still run.
- AMBER writes the metal (`ZNH A 260`) as an `ATOM` record. The *fuzzy* center
  matcher only accepts `HETATM`, so the center is given as a **strict** key
  `ZNH_A260-ZNH_A260` (the dash forces strict mode; the repeated key selects a
  single center).

## Option 1 — CLI (YAML)

From the repository root:

```bash
qp run -c example/amber/config.yaml
```

## Option 2 — Python

```bash
python example/amber/extract_clusters.py
```

## Outputs (`output/small_amber/`)

- `small_amber.pdb` — normalized working copy; `small_amber_original.pdb` — original.
- `A260/` — extracted cluster (PDB + `.xyz`), named after the center's chain+resid.
- `charge.csv`, `count.csv` — per-residue charges and residue counts.

## Generate TeraChem inputs with point charges

After extracting clusters, generate the TeraChem input and MM point-charge
embedding file for each cluster.

CLI:

```bash
qp submit -c example/amber/submit.yaml
```

Python:

```bash
python example/amber/generate_terachem.py
```

Both write, under `output/small_amber/A260/wpbeh/`:

- `qmscript.in` — TeraChem input (includes `pointcharges ptchrges.xyz`)
- `ptchrges.xyz` — MM point charges within the cutoff of the QM centroid
- `jobscript.sh` — scheduler submission script

Notes:

- Point charges are taken from the untouched AMBER original
  (`small_amber_original.pdb`). Its ff14SB names (`HIP`/`HIE`/`HID`) match the
  built-in charge dictionary — the normalized `HIS` copy would not, so Protoss
  is not needed here.
- Metal oxidation state and spin multiplicity come from `ox_info.csv`
  (`total charge = summed residue charges + oxidation`); edit it for your system.
- Non-standard residues/ligands absent from ff14SB (here `NNM`, `SUB`, the `ZNH`
  metal) are excluded from the embedding with a warning. Supply their charges via
  a JSON file and set `charge_embedding_charges` to include them.
