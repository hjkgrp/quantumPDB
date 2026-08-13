Input Formats
=============

QuantumPDB accepts several input formats for specifying which protein structures
to process. All inputs are specified in the YAML configuration file.

Single PDB Code
----------------

Process a single structure by its PDB accession code:

.. code-block:: yaml

   input: 1OS7

The structure is automatically downloaded from the RCSB PDB.

Local PDB File
--------------

Process a local PDB file by providing the path:

.. code-block:: yaml

   input: /path/to/structure.pdb

List of PDB Codes
------------------

Process multiple structures by listing PDB codes:

.. code-block:: yaml

   input: [1OS7, 1FYZ, 1PHM]

CSV Input File
--------------

For high-throughput processing, use a CSV file. This is the recommended format
when processing many structures or when per-PDB parameters are needed.

.. code-block:: yaml

   input: proteins.csv

**Required column:**

- ``pdb_id`` --- PDB accession code or path to a local PDB file

**Optional columns:**

- ``center`` --- Center residue definition for this specific PDB (overrides
  the ``center_residues`` parameter in the YAML config)
- ``force_include_residues`` --- Specific protein residues to force-include for
  this PDB (overrides the ``force_include_residues`` parameter in the YAML
  config). See :ref:`force-include-remove-residues-syntax`.
- ``force_remove_residues`` --- Specific protein residues to force-exclude for
  this PDB (overrides the ``force_remove_residues`` parameter in the YAML
  config). See :ref:`force-include-remove-residues-syntax`.
- ``oxidation`` --- Metal oxidation state (required for ``qp submit``)
- ``multiplicity`` --- Spin multiplicity (required for ``qp submit``)

**Example CSV:**

.. code-block:: text

   pdb_id,center,force_include_residues,force_remove_residues,oxidation,multiplicity
   1OS7,FE,HIS_A123,,3,6
   1FYZ,FE2_A5001-FE2_A5002,,ASN_A45,6,11
   1PHM,CU_A357-CU_A358,HIS_A93-GLU_A45,TYR_A200,4,3

Center Residue Syntax
---------------------

The center residue defines which atoms are at the core of the QM cluster model.
The syntax differs between YAML and CSV formats.

YAML Format (``center_residues``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use brackets in the YAML config for **general mode**, which matches all HETATM
records with the given residue name(s):

.. code-block:: yaml

   # Single residue type
   center_residues: [FE]

   # Multiple residue types (generates separate clusters for each)
   center_residues: [FE, FE2]

General mode only matches HETATM records, not protein ATOM records. This
prevents generating a cluster for every instance of a common amino acid.

CSV Format (``center`` column)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CSV ``center`` column supports **specific mode** with exact residue
identification:

.. code-block:: text

   # Single residue (name only, general mode)
   FE

   # Specific residue (name_chain+number)
   FE_A199

   # Merged centers (dash-separated)
   CU_A357-CU_A358

   # Complex merged center
   FE_A155-OXY_A555-HEM_A155-HIS_A93

The format for specific residues is ``RESNAME_CHAINID``, where ``RESNAME`` is
the three-letter residue code and ``CHAINID`` is the chain letter followed by
the residue number (e.g., ``FE_A199`` for iron at position 199 on chain A).

**Merged centers** treat multiple residues as a single center. This is used for:

- Multi-metallic active sites (e.g., ``CU_A357-CU_A358`` for dicopper)
- Heme groups (e.g., ``FE_A155-HEM_A155-OXY_A555-HIS_A93``)
- Oligomeric substrates where you want the cluster centered on a subset

.. important::

   When both ``center_residues`` (YAML) and a ``center`` column (CSV) are
   provided, the **CSV value takes priority** for each PDB.

Merging Nearby Centers
^^^^^^^^^^^^^^^^^^^^^^

For multi-metallic sites where centers are close together, use the
``merge_distance_cutoff`` parameter to automatically merge centers within a
given distance:

.. code-block:: yaml

   center_residues: [FE2]
   merge_distance_cutoff: 4.0

This merges any ``FE2`` atoms within 4.0 Å of each other into a single center.
Set the cutoff slightly larger than the metal--metal distance. Examples:

- **MMO diiron** (Fe--Fe ~3.4 Å): ``merge_distance_cutoff: 4.0``
- **PHM dicopper** (Cu--Cu ~11 Å): ``merge_distance_cutoff: 12.0``

.. _force-include-remove-residues-syntax:

Force-Include / Force-Remove Residues Syntax
----------------------------------------------

``force_include_residues`` force-includes specific protein residues in the QM
cluster, even if they lie beyond the grown coordination spheres.
``force_remove_residues`` does the opposite: it force-excludes specific
residues, even if sphere growth, ``additional_ligands``, or
``force_include_residues`` would otherwise include them (if the same residue
appears in both, removal wins; the cluster's center residue can't be
removed this way — a matching entry is ignored with a warning instead).
Both share the same syntax. Unlike ``additional_ligands`` (which matches
by residue name), each entry must identify one exact residue with the
same ``RESNAME_CHAINID`` format used for specific-mode centers (e.g.
``HIS_A123``), since a resname like ``HIS`` alone can't distinguish which
occurrence is meant.

YAML Format (``force_include_residues`` / ``force_remove_residues``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use YAML list syntax, one ``RESNAME_CHAINID`` token per residue:

.. code-block:: yaml

   force_include_residues: [HIS_A123, GLU_A45]
   force_remove_residues: [TYR_A200]

Applied identically to every PDB in the run — most useful for a
single-structure run, since chain/residue numbers are structure-specific.

CSV Format (``force_include_residues`` / ``force_remove_residues`` columns)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A CSV cell can't hold a YAML list, so multiple residues are dash-separated
within the cell instead, as shown in the example CSV above (e.g.
``HIS_A93-GLU_A45``). Leave a cell blank for PDBs with nothing to
force-include/exclude.

.. important::

   When both the YAML list and its CSV column are provided for either
   parameter, the **CSV value takes priority** for every PDB in the run,
   exactly like ``center``/``center_residues``. This applies
   independently to ``force_include_residues`` and
   ``force_remove_residues``.

Configuration File Structure
-----------------------------

The YAML configuration file has three sections corresponding to the three CLI
commands. You only need to include the parameters relevant to the command you
are running.

.. code-block:: yaml

   #######
   # RUN #
   #######
   input: proteins.csv
   output_dir: output
   modeller: true
   protoss: true
   coordination: true
   center_residues: [FE]
   number_of_spheres: 2

   ##########
   # SUBMIT #
   ##########
   method: wpbeh
   basis: lacvps_ecp
   create_jobs: true

   ###########
   # ANALYZE #
   ###########
   charge_scheme: Hirshfeld
   calc_charge_schemes: true

See :doc:`configuration` for the complete parameter reference.
