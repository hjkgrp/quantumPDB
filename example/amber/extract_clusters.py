"""Extract QM coordination-sphere clusters from a user-supplied AMBER PDB.

This is the programmatic equivalent of ``qp run -c config.yaml`` (in this same
directory) for an AMBER ff14SB structure. It:

1. Normalizes AMBER protonation-state residue names (HIE/HIP/HID, CYX, ASH/GLH,
   LYN, ...) to canonical PDB names so BioPython recognizes them, while keeping
   an untouched copy of the original for reference/point charges.
2. Skips Modeller/Protoss (the structure is already modeled and protonated).
3. Extracts coordination spheres around the Zn center and writes the cluster
   PDBs/XYZ plus charge.csv and count.csv.

Run from anywhere:
    python example/amber/extract_clusters.py
"""

import os

from qp.structure import setup
from qp.cluster import spheres
from qp.cluster.spheres import CenterResidue

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_PDB = os.path.join(HERE, "pdb", "small_amber.pdb")
OUTPUT_DIR = os.path.join(HERE, "output")

# AMBER writes the Zn as an ATOM record (blank hetflag), which the *fuzzy*
# center matcher rejects. Use a *strict* RESNAME_CHAINRESID key (the dash forces
# strict mode; repeating the key selects a single center).
CENTER = "ZNH_A260-ZNH_A260"


def main():
    pdb = os.path.splitext(os.path.basename(INPUT_PDB))[0]

    # Preserve the original AMBER file and write a canonical-named working copy
    # to output/<pdb>/<pdb>.pdb (the original goes to <pdb>_original.pdb).
    working_pdb, original_pdb = setup.prepare_local_pdb(pdb, INPUT_PDB, OUTPUT_DIR)
    print(f"> Original preserved: {original_pdb}")
    print(f"> Normalized working copy: {working_pdb}")

    struct_out = os.path.join(OUTPUT_DIR, pdb)
    cluster_paths = spheres.extract_clusters(
        working_pdb,
        struct_out,
        CenterResidue(CENTER),
        sphere_count=2,
        first_sphere_radius=4.0,
        max_atom_count=750,
        smooth_method="dummy_atom",
        mean_distance=3,
        capping=1,
        charge=True,
        count=True,
        xyz=True,
    )

    print("> Extracted cluster directories:")
    for path in cluster_paths:
        print(f"    {path}")
    print(f"> Charges:  {os.path.join(struct_out, 'charge.csv')}")
    print(f"> Counts:   {os.path.join(struct_out, 'count.csv')}")


if __name__ == "__main__":
    main()
