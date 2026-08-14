"""Polymer A/C/G/U and MGT charge / O3' capping when absent from Protoss ligands."""

from __future__ import annotations

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from qp.cluster import spheres


def _add_atom(res, name, element, coord, serial):
    atom = Atom(
        name,
        np.asarray(coord, dtype=float),
        0.0,
        1.0,
        " ",
        name.rjust(4)[:4],
        serial,
        element,
    )
    res.add(atom)
    return serial + 1


def _build_rna_a(origin=(0.0, 0.0, 0.0)):
    """Minimal A with terminal phosphate (OP1/OP2 CN=1) and dangling O3'.

    Coordinates keep non-bonded atoms outside the 1.8 A valence search.
    """
    structure = Structure("test")
    model = Model(0)
    chain = Chain("D")
    structure.add(model)
    model.add(chain)

    o = np.asarray(origin, dtype=float)
    res = Residue((" ", 2, " "), "A", " ")
    serial = 1
    serial = _add_atom(res, "P", "P", o + (0.0, 0.0, 0.0), serial)
    serial = _add_atom(res, "OP1", "O", o + (1.55, 0.0, 0.0), serial)
    serial = _add_atom(res, "OP2", "O", o + (-0.78, 1.34, 0.0), serial)
    serial = _add_atom(res, "O5'", "O", o + (-0.78, -1.34, 0.0), serial)
    # Sugar displaced so O3' only sees C3' within 1.8 A
    serial = _add_atom(res, "C5'", "C", o + (-2.2, -1.8, 0.0), serial)
    serial = _add_atom(res, "C4'", "C", o + (-3.5, -1.0, 0.0), serial)
    serial = _add_atom(res, "O4'", "O", o + (-4.5, -1.5, 0.8), serial)
    serial = _add_atom(res, "C1'", "C", o + (-5.2, -0.5, 0.0), serial)
    serial = _add_atom(res, "C3'", "C", o + (-3.5, 0.5, -0.5), serial)
    serial = _add_atom(res, "O3'", "O", o + (-3.5, 1.9, -0.5), serial)
    chain.add(res)
    return structure, res


def _build_mgt(origin=(0.0, 0.0, 0.0)):
    """Minimal MGT: three terminal phosphate pairs + N7, well separated."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("D")
    structure.add(model)
    model.add(chain)

    o = np.asarray(origin, dtype=float)
    res = Residue(("H_MGT", 101, " "), "MGT", " ")
    serial = 1
    for i, tag in enumerate("ABG"):
        p = o + (i * 6.0, 0.0, 0.0)
        serial = _add_atom(res, f"P{tag}", "P", p, serial)
        serial = _add_atom(res, f"O1{tag}", "O", p + (1.55, 0.0, 0.0), serial)
        serial = _add_atom(res, f"O2{tag}", "O", p + (-0.78, 1.34, 0.0), serial)
    # N7 methylated: N7 bonded to C5, C8, CM7
    serial = _add_atom(res, "N7", "N", o + (0.0, 5.0, 0.0), serial)
    serial = _add_atom(res, "C5", "C", o + (1.2, 5.8, 0.0), serial)
    serial = _add_atom(res, "C8", "C", o + (-1.2, 5.8, 0.0), serial)
    serial = _add_atom(res, "CM7", "C", o + (0.0, 3.6, 0.0), serial)
    chain.add(res)
    return structure, res


def test_polymer_a_phosphate_charge_and_o3_cap():
    structure, res = _build_rna_a()
    center = spheres.CenterResidue("exact:A_D2")
    charge = spheres.compute_charge([{res}], structure, {}, center)
    assert charge == [-1]

    spheres.cap_chains(structure[0], {res}, capping=1, ligand_charge={})
    assert spheres.has_res_atom(res, "HO3'")


def test_polymer_a_skipped_when_in_ligand_oligomer():
    structure, res = _build_rna_a()
    center = spheres.CenterResidue("exact:A_D2")
    ligand_charge = {"G_D1 A_D2 U_D3": -3}
    charge = spheres.compute_charge([{res}], structure, ligand_charge, center)
    assert charge == [0]

    spheres.cap_chains(structure[0], {res}, capping=2, ligand_charge=ligand_charge)
    assert not spheres.has_res_atom(res, "HO3'")


def test_mgt_triphosphate_and_n7_charge():
    structure, res = _build_mgt()
    center = spheres.CenterResidue("exact:MGT_D101")
    # -1 * 3 phosphates + 1 for N7
    charge = spheres.compute_charge([{res}], structure, {}, center)
    assert charge == [-2]


def test_mgt_skipped_when_in_ligands():
    structure, res = _build_mgt()
    center = spheres.CenterResidue("exact:MGT_D101")
    ligand_charge = {"MGT_D101": -5}
    charge = spheres.compute_charge([{res}], structure, ligand_charge, center)
    assert charge == [0]


def test_o3_cap_with_ace_nme_flag():
    structure, res = _build_rna_a()
    spheres.cap_chains(structure[0], {res}, capping=2, ligand_charge={})
    assert spheres.has_res_atom(res, "HO3'")
