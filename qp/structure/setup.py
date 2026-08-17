"""Fetch PDB and perform quality checks"""

import os
import csv
import sys
import yaml
import requests

import shutil
from qp.structure.amber_dict import get_amber_to_pdbcanonical_name
from qp.structure.mmcif_to_pdb import OversizedStructureError, convert_mmcif_to_pdb


def read_config(config_file):
    """Read and parse a YAML configuration file.

    Parameters
    ----------
    config_file : str
        Path to the configuration YAML file.

    Returns
    -------
    dict
        Parsed configuration data with all pipeline parameters.
    """
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)


def parse_input(input, output, center_yaml_residues, force_include_yaml_residues=None, force_remove_yaml_residues=None):
    """Parse input sources and determine center residues for each structure.

    Processes the user-provided input (PDB codes, file paths, or CSV) and
    determines the center residue definitions from either the CSV ``center``
    column or the YAML ``center_residues`` parameter. Also resolves the
    per-PDB ``force_include_residues`` and ``force_remove_residues`` lists
    from either their CSV column or their YAML parameter, using the same
    CSV-takes-priority precedence as centers.

    Parameters
    ----------
    input : str or list
        Input specification: PDB code(s), path to PDB/mmCIF file(s), or path to CSV.
    output : str
        Path to the output directory.
    center_yaml_residues : list
        Center residue definitions from the YAML config file.
    force_include_yaml_residues : list, optional
        Additional protein residue keys (``'RESNAME_CHAINID'``, e.g.
        ``'HIS_A123'``) from the YAML config file, applied to every PDB
        when the CSV has no ``force_include_residues`` column (default None,
        treated as []).
    force_remove_yaml_residues : list, optional
        Protein residue keys (``'RESNAME_CHAINID'``) to force-exclude, from
        the YAML config file, applied to every PDB when the CSV has no
        ``force_remove_residues`` column (default None, treated as []).

    Returns
    -------
    tuple of (list, list, list, list)
        ``(pdb_all, center_residues, force_include_residues,
        force_remove_residues)`` where ``pdb_all`` is a list of
        ``(pdb_id, pdb_path, source)`` tuples, ``center_residues`` is a
        list of center definition strings, and ``force_include_residues``/
        ``force_remove_residues`` are each a list (one entry per PDB) of
        lists of ``'RESNAME_CHAINID'`` residue keys to force-include/exclude
        for that PDB. see :func:`get_pdbs` for the
        ``source`` schema.

    Raises
    ------
    SystemExit
        If no center residues are provided in either the CSV or YAML.
    """

    # If the input was a path or a single PDB, convert it to a list
    if isinstance(input, str):
        input = [input]
    output = os.path.abspath(output)
    center_csv_residues = get_centers(input)
    force_include_csv_residues = get_force_include_residues(input)
    force_remove_csv_residues = get_force_remove_residues(input)
    pdb_all = get_pdbs(input, output)

    # Determine how the user has decided to provide the center residues
    if center_csv_residues:
        use_csv_centers = True
        print("> Using residues from the input csv as centers\n")
    elif center_yaml_residues:
        use_csv_centers = False
        print("> Using residues from the input yaml as centers.\n")
    else:
        print("> No center residues were provided in the input csv or the yaml.\n")
        exit()

    # Determine the center residues for the current PDB
    if use_csv_centers:
        center_residues = center_csv_residues
    else:
        if not center_yaml_residues:
            sys.exit("> No more center residues available.")
        center_residues = center_yaml_residues

    # Determine the force-included residues for each PDB
    if force_include_csv_residues:
        print("> Using force-included residues from the input csv\n")
        force_include_residues = [[t for t in row.split("-") if t] for row in force_include_csv_residues]
    else:
        force_include_residues = [list(force_include_yaml_residues or []) for _ in pdb_all]

    # Determine the force-excluded residues for each PDB
    if force_remove_csv_residues:
        print("> Using force-excluded residues from the input csv\n")
        force_remove_residues = [[t for t in row.split("-") if t] for row in force_remove_csv_residues]
    else:
        force_remove_residues = [list(force_remove_yaml_residues or []) for _ in pdb_all]

    return pdb_all, center_residues, force_include_residues, force_remove_residues


def fetch_pdb(pdb, out):
    """
    Fetch a structure for a PDB code, falling back to mmCIF when needed.

    Tries the classic PDB download first. If RCSB returns 404 (common for
    entries with ``pdb_format_compatible = N``), downloads the mmCIF and
    converts it to ``out``.

    Parameters
    ----------
    pdb: str
        PDB code
    out: str
        Path to output PDB file

    Returns
    -------
    str
        ``"pdb"`` if a classic PDB was downloaded, or ``"mmcif"`` if the
        file was obtained via mmCIF conversion.

    Raises
    ------
    ValueError
        If neither PDB nor mmCIF is available for the ID.
    IOError
        If there is a network or server-side error.
    """
    pdb = pdb.strip()
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    pdb_url = f"https://files.rcsb.org/download/{pdb}.pdb"
    try:
        r = requests.get(pdb_url, timeout=30)
    except requests.exceptions.RequestException as e:
        raise IOError(f"Could not connect to the server. Details: {e}")

    if r.status_code == 200:
        with open(out, "w") as f:
            f.write(r.text)
        return "pdb"
    if r.status_code != 404:
        raise IOError(f"Server returned an error with status code {r.status_code}.")

    # Classic PDB unavailable — try mmCIF (entries marked pdb_format_compatible=N)
    print(f"> PDB unavailable for {pdb}; fetching mmCIF")
    cif_url = f"https://files.rcsb.org/download/{pdb}.cif"
    try:
        r_cif = requests.get(cif_url, timeout=60)
    except requests.exceptions.RequestException as e:
        raise IOError(f"Could not connect to the server. Details: {e}")

    if r_cif.status_code == 404:
        raise ValueError(f"PDB ID '{pdb}' is not valid or does not exist.")
    if r_cif.status_code != 200:
        raise IOError(f"Server returned an error with status code {r_cif.status_code}.")

    cif_path = os.path.join(os.path.dirname(out), f"{pdb}.cif")
    with open(cif_path, "w") as f:
        f.write(r_cif.text)
    print(f"> Converting mmCIF → PDB ({cif_path})")
    convert_mmcif_to_pdb(cif_path, out)
    return "mmcif"


def ensure_structure_pdb(pdb_id, pdb_path, source_cif=None):
    """Ensure ``pdb_path`` exists as a classic PDB, converting mmCIF if needed.

    Parameters
    ----------
    pdb_id : str
        Structure identifier (PDB code or local basename).
    pdb_path : str
        Target classic PDB path used by the rest of the pipeline.
    source_cif : str or None
        Optional local mmCIF path to convert when ``pdb_path`` is missing.

    Returns
    -------
    str
        One of ``"exists"``, ``"converted"``, ``"pdb"``, or ``"mmcif"``.
    """
    pdb_path = os.path.abspath(pdb_path)
    if os.path.isfile(pdb_path):
        return "exists"

    if source_cif:
        print(f"> Converting mmCIF → PDB ({source_cif})")
        convert_mmcif_to_pdb(source_cif, pdb_path)
        return "converted"

    return fetch_pdb(pdb_id, pdb_path)

def convert_amber_to_pdb(input_pdb, output_pdb):
    """Rewrite AMBER ff14SB residue names to canonical PDB residue names.

    Parameters
    ----------
    input_pdb: str
        Path to the input PDB file
    output_pdb: str
        Path to the output PDB file

    Returns
    -------
    changed: set
        Set of tuples (resname, chain, resid) that were changed
    """
    amber_to_pdbcanonical = get_amber_to_pdbcanonical_name()
    changed = set()
    with open(input_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                if resname in amber_to_pdbcanonical:
                    line = line[:17] + amber_to_pdbcanonical[resname] + line[20:]
                    changed.add((resname,chain,resid))
            outfile.write(line)
    return changed

def prepare_local_pdb(pdb, user_pdb_path, output_path):
    
    struct_dir = os.path.join(output_path, pdb)
    os.makedirs(struct_dir, exist_ok=True)

    original = os.path.join(struct_dir, f"{pdb}_original.pdb")
    target = os.path.join(struct_dir, f"{pdb}.pdb")
    
    shutil.copy(user_pdb_path, original)
    changed = convert_amber_to_pdb(original, target)
    if changed:
        summary = ", ".join(f"{resname} {chain}{resid}" for resname, chain, resid in sorted(changed))
        print(f"> Renamed AMBER residues in {pdb}: {summary}")
    
    return target, original

def get_pdbs(input_path, output_path):
    """
    Parses the input PDBs and returns a list of tuples

    Parameters
    ----------
    input_path: list
        List of PDB codes, paths to PDB/mmCIF files, or path to CSV file
    output_path: str
        Path to output directory

    Returns
    -------
    pdb_all
        list of tuple
    One ``(pdb_id, pdb_path, source)`` tuple per structure. ``pdb_path`` is
    the working PDB the rest of the pipeline consumes. ``source`` records
    is one of:
    * ``None`` — fetched by PDB code, or entry from a CSV/txt batch.
    * ``{"type": "cif", "path": <abs mmCIF path>}`` — local mmCIF that is
      converted to ``pdb_path`` during the run.
    * ``{"type": "amber", "path": <abs original PDB path>}`` — local PDB
      whose AMBER residue names were normalized into ``pdb_path``; the
      untouched original is preserved at ``path`` for point-charge
      embedding.

    Notes
    -----
    Store input PDBs as a tuple of
    Parsed ID (PDB code or filename)
    Path to PDB file (existing or to download/convert)
    Optional local mmCIF source path
    """
    pdb_all = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            ext_lower = ext.lower()
            if ext_lower == ".pdb":
                target, original = prepare_local_pdb(pdb, pdb_id, output_path)
                pdb_all.append((pdb, target, {"type": "amber", "path": original}))
            elif ext_lower in (".cif", ".mmcif"):
                target = os.path.join(output_path, pdb, f"{pdb}.pdb")
                pdb_all.append((pdb, target, {"type": "cif", "path": os.path.abspath(pdb_id)}))
            elif ext_lower == ".csv":
                with open(pdb_id, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        pdb = row['pdb_id']
                        pdb_all.append((pdb, os.path.join(output_path, pdb, f"{pdb}.pdb"), None))
            else:
                with open(pdb_id, "r") as f:
                    pdb_all.extend(
                        [
                            (pdb, os.path.join(output_path, pdb, f"{pdb}.pdb"), None)
                            for pdb in f.read().splitlines()
                            if pdb.strip()
                        ]
                    )
        else:
            pdb_all.append((pdb_id, os.path.join(output_path, pdb_id, f"{pdb_id}.pdb"), None))

    return pdb_all


def get_centers(input_path):
    """
    Parses the input centers and returns a list of the center residue PDB IDs.

    Parameters
    ----------
    input_path: list
        List of pdbs or the input csv file.

    Returns
    -------
    centers: list
        List PDB IDs for user-defined center residues.
    """
    
    centers = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            if ext == ".csv":
                input_csv = pdb_id
                with open(input_csv, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    # Check if 'center' column exists
                    if 'center' not in reader.fieldnames:
                        print(f"> WARNING: The 'center' option is not being used in {input_csv}. Returning empty list.")
                        return []
                    # If the 'center' column exists, proceed to collect centers
                    for row in reader:
                        center = row.get('center', None)
                        if center:
                            centers.append(center)
    return centers


def get_force_include_residues(input_path):
    """
    Parses the input csv's ``force_include_residues`` column, one entry per row.

    Unlike ``get_centers``, blank cells append an empty string rather than
    being skipped, since ``force_include_residues`` is optional per PDB and the
    returned list must stay aligned with the CSV's rows.

    Parameters
    ----------
    input_path: list
        List of pdbs or the input csv file.

    Returns
    -------
    force_include_residues: list
        List of raw ``force_include_residues`` cell values (dash-separated
        ``'RESNAME_CHAINID'`` tokens, or ``""``), one per CSV row. Empty
        list if no CSV with an ``force_include_residues`` column is present.
    """

    force_include_residues = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            if ext == ".csv":
                input_csv = pdb_id
                with open(input_csv, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    # Check if 'force_include_residues' column exists
                    if 'force_include_residues' not in reader.fieldnames:
                        print(f"> WARNING: The 'force_include_residues' option is not being used in {input_csv}. Returning empty list.")
                        return []
                    # If the column exists, proceed to collect its values
                    for row in reader:
                        force_include_residues.append(row.get('force_include_residues', None) or "")
    return force_include_residues


def get_force_remove_residues(input_path):
    """
    Parses the input csv's ``force_remove_residues`` column, one entry per row.

    Unlike ``get_centers``, blank cells append an empty string rather than
    being skipped, since ``force_remove_residues`` is optional per PDB and the
    returned list must stay aligned with the CSV's rows.

    Parameters
    ----------
    input_path: list
        List of pdbs or the input csv file.

    Returns
    -------
    force_remove_residues: list
        List of raw ``force_remove_residues`` cell values (dash-separated
        ``'RESNAME_CHAINID'`` tokens, or ``""``), one per CSV row. Empty
        list if no CSV with an ``force_remove_residues`` column is present.
    """

    force_remove_residues = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            if ext == ".csv":
                input_csv = pdb_id
                with open(input_csv, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    # Check if 'force_remove_residues' column exists
                    if 'force_remove_residues' not in reader.fieldnames:
                        print(f"> WARNING: The 'force_remove_residues' option is not being used in {input_csv}. Returning empty list.")
                        return []
                    # If the column exists, proceed to collect its values
                    for row in reader:
                        force_remove_residues.append(row.get('force_remove_residues', None) or "")
    return force_remove_residues
