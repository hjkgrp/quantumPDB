"""Tests for AMBER ff14SB residue-name normalization of user-supplied PDBs."""

import os

from Bio.PDB.Polypeptide import is_aa

from qp.structure import setup
from qp.structure.amber_dict import (
    get_amber_to_pdbcanonical_name,
    get_resname_dict,
)
from qp.cluster import spheres
from qp.cluster.spheres import CenterResidue

# A real AMBER-prepared structure (contains HIE/HIP among standard residues,
# plus custom non-AMBER names like ZNH/SUB that must be left untouched).
SMALL_AMBER = os.path.join(os.path.dirname(__file__), "data", "small_amber.pdb")


def _residue_atoms(pdb_path):
    """Map ``(chain, resid) -> [resname, {atom_names}]`` via fixed-column parsing.

    Line-based (not Bio.PDB) to stay robust to the fixture's PQR-style
    charge/radius columns.
    """
    residues = {}
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                atom = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                resid = line[22:26].strip()
                key = (chain, resid)
                if key not in residues:
                    residues[key] = [resname, set()]
                residues[key][1].add(atom)
    return residues


def _atom(serial, name, resname, chain, resseq, x, y, z,
          occ=1.00, b=20.00, element="C", record="ATOM"):
    """Build a column-accurate PDB ATOM/HETATM line (cols per PDB v3.3)."""
    return (
        f"{record:<6}"        # 1-6   record name
        f"{serial:>5}"        # 7-11  serial
        f" "                  # 12    space
        f"{name:<4}"          # 13-16 atom name
        f" "                  # 17    altLoc
        f"{resname:>3}"       # 18-20 resName
        f" "                  # 21    space
        f"{chain:1}"          # 22    chainID
        f"{resseq:>4}"        # 23-26 resSeq
        f" "                  # 27    iCode
        f"   "                # 28-30 spaces
        f"{x:8.3f}{y:8.3f}{z:8.3f}"  # 31-54 coords
        f"{occ:6.2f}{b:6.2f}"        # 55-66 occ / B-factor
        f"          "         # 67-76 spaces
        f"{element:>2}"       # 77-78 element
        "\n"
    )


# A small structure mixing AMBER protonation-state names, a standard residue,
# a non-AMBER HETATM, and non-ATOM records that must survive conversion.
_AMBER_PDB = (
    "HEADER    TEST AMBER STRUCTURE\n"
    + _atom(1, "N",  "HIP", "A", 10, 11.104, 13.207, 10.000)
    + _atom(2, "SG", "CYX", "A", 11, 12.000, 14.000, 10.000, element="S")
    + _atom(3, "OE1", "GLH", "A", 12, 13.000, 15.000, 10.000, element="O")
    + _atom(4, "CA", "ALA", "A", 13, 14.000, 16.000, 10.000)
    + _atom(5, "FE", "FE",  "A", 99, 15.000, 17.000, 10.000, element="FE",
            record="HETATM")
    + "TER\n"
    + "END\n"
)

_NO_AMBER_PDB = (
    "HEADER    PLAIN STRUCTURE\n"
    + _atom(1, "N",  "HIS", "A", 10, 11.104, 13.207, 10.000)
    + _atom(2, "CA", "ALA", "A", 11, 12.000, 14.000, 10.000)
    + "END\n"
)


# --------------------------------------------------------------------------- #
# Layer 1a: the conversion dictionary
# --------------------------------------------------------------------------- #

def test_amber_to_canonical_mapping():
    m = get_amber_to_pdbcanonical_name()
    assert m["HIP"] == "HIS"
    assert m["HID"] == "HIS"
    assert m["HIE"] == "HIS"
    assert m["CYX"] == "CYS"
    assert m["CYM"] == "CYS"
    assert m["ASH"] == "ASP"
    assert m["GLH"] == "GLU"
    assert m["LYN"] == "LYS"


def test_amber_to_canonical_is_variant_keyed():
    # It maps AMBER variant -> canonical, not canonical -> itself.
    m = get_amber_to_pdbcanonical_name()
    for canonical in get_resname_dict():
        assert canonical not in m
    # Every variant listed in the grouped dict is present in the flat map.
    for canonical, variants in get_resname_dict().items():
        for variant in variants:
            assert m[variant] == canonical


# --------------------------------------------------------------------------- #
# Layer 1b: convert_amber_to_pdb
# --------------------------------------------------------------------------- #

def test_convert_renames_amber_residues(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(_AMBER_PDB)
    out = tmp_path / "out.pdb"

    setup.convert_amber_to_pdb(str(src), str(out))
    text = out.read_text()

    # AMBER names replaced with canonical names in the resName column.
    assert "HIS A" in text
    assert "CYS A" in text
    assert "GLU A" in text
    # Untouched residues stay as-is.
    assert "ALA A" in text
    # Non-AMBER HETATM residue name is left alone.
    assert "FE " in text


def test_convert_preserves_non_atom_records(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(_AMBER_PDB)
    out = tmp_path / "out.pdb"

    setup.convert_amber_to_pdb(str(src), str(out))
    text = out.read_text()

    # HEADER / TER / END must survive (guards against dropping non-ATOM lines).
    assert text.startswith("HEADER    TEST AMBER STRUCTURE")
    assert "TER\n" in text
    assert text.rstrip().endswith("END")
    # HETATM line preserved as a HETATM record.
    assert any(line.startswith("HETATM") for line in text.splitlines())


def test_convert_preserves_columns(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(_AMBER_PDB)
    out = tmp_path / "out.pdb"

    setup.convert_amber_to_pdb(str(src), str(out))

    for line in out.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            assert len(line) >= 66
            # Occupancy and B-factor columns remain independently parseable.
            assert abs(float(line[54:60]) - 1.00) < 1e-6
            assert abs(float(line[60:66]) - 20.00) < 1e-6


def test_convert_returns_changed_residues(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(_AMBER_PDB)
    out = tmp_path / "out.pdb"

    changed = setup.convert_amber_to_pdb(str(src), str(out))

    changed_names = {resname for resname, _chain, _resid in changed}
    assert changed_names == {"HIP", "CYX", "GLH"}


def test_convert_no_amber_leaves_content_and_empty_set(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(_NO_AMBER_PDB)
    out = tmp_path / "out.pdb"

    changed = setup.convert_amber_to_pdb(str(src), str(out))

    assert changed == set()
    assert out.read_text() == _NO_AMBER_PDB


# --------------------------------------------------------------------------- #
# Layer 1c: prepare_local_pdb
# --------------------------------------------------------------------------- #

def test_prepare_local_pdb_preserves_original_and_normalizes(tmp_path):
    user_pdb = tmp_path / "myprot.pdb"
    user_pdb.write_text(_AMBER_PDB)
    output = tmp_path / "out"

    target, original = setup.prepare_local_pdb("myprot", str(user_pdb), str(output))

    # Both artifacts land under {output}/{pdb}/ with the expected names.
    assert target == os.path.join(str(output), "myprot", "myprot.pdb")
    assert original == os.path.join(str(output), "myprot", "myprot_original.pdb")
    assert os.path.isfile(target)
    assert os.path.isfile(original)

    # Original is byte-for-byte the user's input (AMBER names retained).
    assert open(original).read() == _AMBER_PDB
    assert "HIP A" in open(original).read()

    # Working copy is normalized to canonical names.
    working = open(target).read()
    assert "HIS A" in working
    assert "HIP" not in working


# --------------------------------------------------------------------------- #
# Layer 2: get_pdbs provenance
# --------------------------------------------------------------------------- #

def test_get_pdbs_local_pdb_tagged_amber(tmp_path):
    user_pdb = tmp_path / "1abc.pdb"
    user_pdb.write_text(_AMBER_PDB)
    output = tmp_path / "out"

    pdb_all = setup.get_pdbs([str(user_pdb)], str(output))

    assert len(pdb_all) == 1
    pdb_id, path, source = pdb_all[0]
    assert pdb_id == "1abc"
    assert path == os.path.join(str(output), "1abc", "1abc.pdb")
    assert source["type"] == "amber"
    assert source["path"] == os.path.join(str(output), "1abc", "1abc_original.pdb")
    # Original preserved and untouched; working copy normalized.
    assert os.path.isfile(source["path"])
    assert "HIP A" in open(source["path"]).read()
    assert "HIS A" in open(path).read()


def test_get_pdbs_fetched_code_has_no_source(tmp_path):
    output = tmp_path / "out"
    pdb_all = setup.get_pdbs(["1lm6"], str(output))

    assert len(pdb_all) == 1
    pdb_id, path, source = pdb_all[0]
    assert pdb_id == "1lm6"
    assert path == os.path.join(str(output), "1lm6", "1lm6.pdb")
    assert source is None


# --------------------------------------------------------------------------- #
# Layer 2b: real AMBER fixture (qp/tests/data/small_amber.pdb)
# --------------------------------------------------------------------------- #

def test_real_amber_motivation_is_aa():
    # Documents why normalization is needed: AMBER names are not standard AAs.
    # HIE/CYX are unrecognized outright; HIP passes loose is_aa but fails the
    # standard=True check that qp.protonate.fix relies on.
    assert not is_aa("HIE")
    assert not is_aa("CYX")
    assert not is_aa("HIP", standard=True)
    # After normalization everything is canonical HIS -> recognized either way.
    assert is_aa("HIS")
    assert is_aa("HIS", standard=True)


def test_real_amber_preserves_all_records(tmp_path):
    out = tmp_path / "converted.pdb"
    setup.convert_amber_to_pdb(SMALL_AMBER, str(out))

    src_lines = open(SMALL_AMBER).read().splitlines()
    out_lines = out.read_text().splitlines()

    # No records dropped or added; atom-line count unchanged.
    assert len(out_lines) == len(src_lines)
    assert sum(l.startswith(("ATOM", "HETATM")) for l in out_lines) == \
           sum(l.startswith(("ATOM", "HETATM")) for l in src_lines)
    # Non-ATOM framing records survive.
    assert any(l.startswith("CRYST1") for l in out_lines)
    assert any(l.startswith("TER") for l in out_lines)
    assert out_lines[-1].startswith("END")


def test_real_amber_no_dict_names_remain(tmp_path):
    out = tmp_path / "converted.pdb"
    setup.convert_amber_to_pdb(SMALL_AMBER, str(out))

    amber_names = set(get_amber_to_pdbcanonical_name())
    src_names = {r[0] for r in _residue_atoms(SMALL_AMBER).values()}
    out_names = {r[0] for r in _residue_atoms(str(out)).values()}

    # Fixture actually exercises the mapping.
    assert amber_names & src_names, "fixture has no AMBER names to convert"
    # After conversion, no AMBER protonation-state name remains.
    assert not (amber_names & out_names)


def test_real_amber_changed_reports_present_variants(tmp_path):
    out = tmp_path / "converted.pdb"
    changed = setup.convert_amber_to_pdb(SMALL_AMBER, str(out))

    changed_names = {resname for resname, _chain, _resid in changed}
    amber_names = set(get_amber_to_pdbcanonical_name())
    src_names = {r[0] for r in _residue_atoms(SMALL_AMBER).values()}

    # Reports exactly the AMBER names that were present in the input.
    assert changed_names == (amber_names & src_names)


def test_real_amber_untouched_nonstandard_residues(tmp_path):
    out = tmp_path / "converted.pdb"
    setup.convert_amber_to_pdb(SMALL_AMBER, str(out))

    out_names = {r[0] for r in _residue_atoms(str(out)).values()}
    # Custom AMBER names outside the dictionary are preserved verbatim.
    for name in ("ZNH", "SUB", "ACE", "NME"):
        assert name in out_names


def test_real_amber_protonation_hydrogens_preserved(tmp_path):
    """The crux for charge counting: renaming keeps the AMBER hydrogens."""
    out = tmp_path / "converted.pdb"
    setup.convert_amber_to_pdb(SMALL_AMBER, str(out))

    before = _residue_atoms(SMALL_AMBER)
    after = _residue_atoms(str(out))

    hip_keys = [k for k, (name, _) in before.items() if name == "HIP"]
    hie_keys = [k for k, (name, _) in before.items() if name == "HIE"]
    assert hip_keys and hie_keys  # fixture sanity

    # HIP (doubly protonated) -> HIS, keeps BOTH HD1 and HE2 => counts as +1.
    for key in hip_keys:
        resname, atoms = after[key]
        assert resname == "HIS"
        assert "HD1" in atoms and "HE2" in atoms

    # HIE (epsilon only) -> HIS, keeps HE2 but not HD1 => neutral.
    for key in hie_keys:
        resname, atoms = after[key]
        assert resname == "HIS"
        assert "HE2" in atoms and "HD1" not in atoms


# --------------------------------------------------------------------------- #
# Layer 3: end-to-end cluster extraction from the normalized AMBER structure
# --------------------------------------------------------------------------- #

def test_extract_clusters_from_amber(tmp_path):
    """Normalize the AMBER structure, then actually extract coordination spheres.

    A strict center key (``ZNH_A260``) is used because AMBER writes the Zn as an
    ATOM record (blank hetflag), which the fuzzy center matcher rejects.
    """
    work = tmp_path / "work.pdb"
    setup.convert_amber_to_pdb(SMALL_AMBER, str(work))

    out = tmp_path / "out"
    out.mkdir()

    center = CenterResidue("ZNH_A260-ZNH_A260")
    cluster_paths = spheres.extract_clusters(
        str(work), str(out), center,
        sphere_count=2, smooth_method="dummy_atom", mean_distance=3,
        charge=True, count=True, capping=1,
    )

    center_dir = os.path.join(str(out), "A260")
    assert cluster_paths == [center_dir]

    # Center sphere + 2 coordination spheres, all non-empty.
    sphere_resnames = set()
    for i in range(3):
        sphere_pdb = os.path.join(center_dir, f"{i}.pdb")
        assert os.path.isfile(sphere_pdb)
        atoms = [l for l in open(sphere_pdb) if l.startswith(("ATOM", "HETATM"))]
        assert atoms, f"sphere {i} is empty"
        sphere_resnames.update(l[17:20].strip() for l in atoms)

    # Zn center captured; renamed His residues flow through as canonical HIS
    # (no AMBER names leak) -- proving normalization reaches clustering.
    assert "ZNH" in sphere_resnames
    assert "HIS" in sphere_resnames
    assert not ({"HIE", "HID", "HIP"} & sphere_resnames)

    # Charge/residue-count bookkeeping was produced.
    assert os.path.isfile(os.path.join(str(out), "charge.csv"))
    assert os.path.isfile(os.path.join(str(out), "count.csv"))
