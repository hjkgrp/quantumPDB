"""Generate TeraChem input files (with MM point charges) for the AMBER example.

Programmatic equivalent of ``qp submit -c submit.yaml``. For every extracted
cluster it writes, under ``output/<pdb>/<center>/<method>/``:

* ``qmscript.in``  - TeraChem input (references ``ptchrges.xyz``)
* ``ptchrges.xyz`` - MM point charges from ff14SB
* ``jobscript.sh`` - scheduler submission script

Point charges are read from the untouched AMBER original
(``<pdb>_original.pdb``), whose ff14SB protonation-state names (HIP/HIE/HID)
match the built-in charge dictionary; Protoss is not required.

If the clusters have not been extracted yet, this script runs the extraction
first (via ``extract_clusters.py`` in this directory).

Run from anywhere:
    python example/amber/generate_terachem.py
"""

import os
import sys

from qp.manager import create

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(HERE, "output")
# CSV providing metal oxidation state + spin multiplicity per structure.
# total_charge = (residue charges from charge.csv) + oxidation.
ELECTRONIC_CSV = os.path.join(HERE, "ox_info.csv")

METHOD = "wpbeh"


def ensure_clusters():
    """Extract clusters first if the expected output is missing."""
    if os.path.isdir(os.path.join(OUTPUT_DIR, "small_amber")):
        return
    print("> Clusters not found; running extraction first")
    sys.path.insert(0, HERE)
    import extract_clusters

    extract_clusters.main()


def main():
    ensure_clusters()

    create.create_jobs(
        ELECTRONIC_CSV,
        OUTPUT_DIR,
        optimization=False,
        basis="lacvps_ecp",
        method=METHOD,
        guess="generate",
        use_charge_embedding=True,
        charge_embedding_cutoff=20,
        charge_embedding_charges=None,   # None = built-in AMBER ff14SB charges
        gpus=1,
        memory="8G",
        scheduler="slurm",
        pcm_radii_file="pcm_radii",
        dielectric=10,
        use_implicit_solvent=False,      # point-charge embedding instead of PCM
    )

    print("\n> Wrote TeraChem inputs under:")
    for root, _dirs, files in sorted(os.walk(OUTPUT_DIR)):
        if root.endswith(os.sep + METHOD):
            for name in sorted(files):
                print(f"    {os.path.join(root, name)}")


if __name__ == "__main__":
    main()
