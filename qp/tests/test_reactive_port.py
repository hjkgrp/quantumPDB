"""Tests for reactive-branch features ported onto hjkgrp main."""

import os

import numpy as np
import pytest
from Bio.PDB import PDBParser, Polypeptide

from qp.cluster.spheres import (
    CenterResidue,
    find_RGP_atoms,
    get_center_residues,
    get_next_neighbors,
    make_res_key,
    voronoi,
)
from qp.manager.create import ligand_in_spheres
from qp.protonate.ligand_prop import collect_RGP_atoms, parse_ligand_name, read_ligands

MISSING_LICENSE = False
try:
    import modeller  # noqa: F401
except Exception:
    MISSING_LICENSE = True

if not MISSING_LICENSE:
    from qp.structure import missing as missing_loops


SAMPLES = os.path.join(os.path.dirname(__file__), "samples")
TAN_SDF = os.path.join(SAMPLES, "3x20", "Protoss", "3x20_ligands.sdf")


def test_collect_rgp_atoms_from_3x20():
    rgp = collect_RGP_atoms(TAN_SDF)
    assert "TAN_A302" in rgp
    assert len(rgp["TAN_A302"]) == 1
    atom_info = next(iter(rgp["TAN_A302"].values()))
    assert np.allclose(atom_info["coord"], [105.6540, 26.7590, 17.5180], atol=1e-4)
    # Bond 5-17 links C5 to R#17 in the TAN block.
    assert np.allclose(
        atom_info["linking_atom_coord"], [105.8300, 25.2830, 17.1380], atol=1e-4
    )


def test_parse_ligand_name_and_read_ligands():
    blocks = read_ligands(TAN_SDF)
    tan_block = next(block for block in blocks if block and block[0].startswith("TAN_"))
    name_id, name = parse_ligand_name(tan_block)
    assert name == "TAN_A302"
    assert name_id == [("TAN", "A", "302")]


def test_find_rgp_atoms_matches_pdb_coords(tmpdir):
    # Minimal PDB with the TAN RGP and linking-atom coordinates from the SDF.
    pdb_text = """\
ATOM      1  C5  TAN A 302     105.830  25.283  17.138  1.00  0.00           C
ATOM      2  R1  TAN A 302     105.654  26.759  17.518  1.00  0.00           R
END
"""
    pdb_path = os.path.join(tmpdir, "tan.pdb")
    with open(pdb_path, "w") as f:
        f.write(pdb_text)
    structure = PDBParser(QUIET=True).get_structure("TAN", pdb_path)
    rgp = collect_RGP_atoms(TAN_SDF)
    find_RGP_atoms(structure, rgp)
    info = next(iter(rgp["TAN_A302"].values()))
    assert info["atom"].get_name() == "R1"
    assert info["linking_atom"].get_name() == "C5"
    assert make_res_key(info["atom"].get_parent()) == "TAN_A302"


def test_include_ligands_mode_3_excludes_nonwater_ligands(tmpdir):
    pdb_path = os.path.join(SAMPLES, "1lm6", "Protoss", "1lm6_protoss.pdb")
    structure = PDBParser(QUIET=True).get_structure("PDB", pdb_path)
    model = structure[0]

    center = CenterResidue("FE")
    neighbors = voronoi(model, center, [], "dummy_atom", str(tmpdir), mean_distance=3)
    centers = get_center_residues(model, center, 0.0)
    assert centers
    _, _residues, spheres = get_next_neighbors(
        centers[0],
        neighbors,
        sphere_count=2,
        ligands=[],
        first_sphere_radius=4.0,
        smooth_method="dummy_atom",
        include_ligands=3,
        mean_distance=3,
    )
    outer = set().union(*spheres[1:]) if len(spheres) > 1 else set()
    for res in outer:
        if not Polypeptide.is_aa(res):
            assert res.get_resname() == "HOH"


@pytest.mark.skipif(MISSING_LICENSE, reason="Modeller license not found")
def test_get_residues_ter_starts_new_chain_bucket(tmpdir):
    pdb_text = """\
REMARK 465     MISSING RESIDUES
REMARK 465   M RES C SSSEQI
REMARK 465     ALA A    3
ATOM      1  N   ALA A   1      11.104  13.207   1.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.000  13.000   1.000  1.00  0.00           C
ATOM      3  C   ALA A   1      13.000  13.000   1.000  1.00  0.00           C
ATOM      4  O   ALA A   1      13.500  14.000   1.000  1.00  0.00           O
ATOM      5  N   ALA A   2      14.000  13.000   1.000  1.00  0.00           N
ATOM      6  CA  ALA A   2      15.000  13.000   1.000  1.00  0.00           C
ATOM      7  C   ALA A   2      16.000  13.000   1.000  1.00  0.00           C
ATOM      8  O   ALA A   2      16.500  14.000   1.000  1.00  0.00           O
TER       9      ALA A   2
ATOM     10  N   ALA B   1      21.104  13.207   1.000  1.00  0.00           N
ATOM     11  CA  ALA B   1      22.000  13.000   1.000  1.00  0.00           C
ATOM     12  C   ALA B   1      23.000  13.000   1.000  1.00  0.00           C
ATOM     13  O   ALA B   1      23.500  14.000   1.000  1.00  0.00           O
END
"""
    pdb_path = os.path.join(tmpdir, "ter.pdb")
    with open(pdb_path, "w") as f:
        f.write(pdb_text)

    AA = missing_loops.define_residues()
    residues = missing_loops.get_residues(pdb_path, AA)
    residues = missing_loops.clean_termini(residues)
    assert len(residues) == 2
    assert [res[0][0] for res in residues[0]] == [1, 2, 3]
    assert residues[0][-1][2] == "R"  # missing ALA A3 appended before TER cut
    assert [res[0][0] for res in residues[1]] == [1]


def test_ligand_in_spheres_reports_missing_monomer(tmpdir, capsys):
    cluster = os.path.join(tmpdir, "cluster")
    os.makedirs(cluster)
    # Only the first monomer is present in sphere 0.
    with open(os.path.join(cluster, "0.pdb"), "w") as f:
        f.write(
            "HETATM    1  C1  GAL A   1      0.000   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
    with open(os.path.join(cluster, "1.pdb"), "w") as f:
        f.write("END\n")

    assert ligand_in_spheres("GAL_A1 GAL_A2", cluster, 1) is False
    captured = capsys.readouterr().out
    assert "GAL_A2 not found" in captured
    assert "GAL_A1 GAL_A2" in captured
