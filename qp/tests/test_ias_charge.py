"""Sphere-0 hetero / IAS formal charges when Protoss omits them from ligands.sdf."""

from __future__ import annotations

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.NeighborSearch import NeighborSearch
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


def _build_ias_isopeptide():
    """IAS with OXT carboxylate and CG isopeptide-linked to HIS N (no OD2)."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("C")
    structure.add(model)
    model.add(chain)

    # Preceding residue so IAS is not treated as the chain N-terminus.
    # Keep PRO CA outside the 1.8 A N search so only the peptide C counts.
    pro = Residue((" ", 3, " "), "PRO", " ")
    serial = 1
    serial = _add_atom(pro, "N", "N", (-4.5, 1.0, 0.0), serial)
    serial = _add_atom(pro, "CA", "C", (-3.2, 0.8, 0.0), serial)
    serial = _add_atom(pro, "C", "C", (-1.35, 0.0, 0.0), serial)
    chain.add(pro)

    ias = Residue(("H_IAS", 4, " "), "IAS", " ")
    serial = _add_atom(ias, "N", "N", (0.0, 0.0, 0.0), serial)
    serial = _add_atom(ias, "H2", "H", (0.0, 1.0, 0.0), serial)
    serial = _add_atom(ias, "CA", "C", (1.5, 0.0, 0.0), serial)
    serial = _add_atom(ias, "C", "C", (2.5, 1.2, 0.0), serial)
    serial = _add_atom(ias, "O", "O", (2.5, 2.4, 0.0), serial)
    serial = _add_atom(ias, "OXT", "O", (3.7, 0.6, 0.0), serial)
    serial = _add_atom(ias, "CB", "C", (1.5, -1.5, 0.0), serial)
    serial = _add_atom(ias, "CG", "C", (1.5, -3.0, 0.0), serial)
    serial = _add_atom(ias, "OD1", "O", (2.7, -3.6, 0.0), serial)
    chain.add(ias)

    his = Residue((" ", 5, " "), "HIS", " ")
    serial = _add_atom(his, "N", "N", (0.3, -3.6, 0.0), serial)
    serial = _add_atom(his, "CA", "C", (-1.0, -4.2, 0.0), serial)
    chain.add(his)
    return structure, ias, his


def _build_ias_free_sidechain():
    """IAS with free OD1/OD2 carboxylate (no isopeptide) and OXT."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("C")
    structure.add(model)
    model.add(chain)

    ias = Residue(("H_IAS", 4, " "), "IAS", " ")
    serial = 1
    serial = _add_atom(ias, "N", "N", (0.0, 0.0, 0.0), serial)
    serial = _add_atom(ias, "CA", "C", (1.5, 0.0, 0.0), serial)
    serial = _add_atom(ias, "C", "C", (2.5, 1.2, 0.0), serial)
    serial = _add_atom(ias, "O", "O", (2.5, 2.4, 0.0), serial)
    serial = _add_atom(ias, "OXT", "O", (3.7, 0.6, 0.0), serial)
    serial = _add_atom(ias, "CB", "C", (1.5, -1.5, 0.0), serial)
    serial = _add_atom(ias, "CG", "C", (1.5, -3.0, 0.0), serial)
    serial = _add_atom(ias, "OD1", "O", (2.7, -3.6, 0.0), serial)
    serial = _add_atom(ias, "OD2", "O", (0.3, -3.6, 0.0), serial)
    chain.add(ias)
    return structure, ias


def test_ias_isopeptide_gets_oxt_minus_one_not_sidechain():
    structure, ias, his = _build_ias_isopeptide()
    tree = NeighborSearch(list(structure.get_atoms()))
    assert spheres.ias_cg_isopeptide_linked(ias, tree)
    assert spheres.hetero_residue_formal_charge(ias, tree, set()) == -1

    # Fuzzy center: sphere-0 heteros land in ligand_charge via compute_charge.
    center = spheres.CenterResidue("ADN_IAS")
    ligand_charge = {}
    spheres.compute_charge([{ias}], structure, ligand_charge, center)
    assert ligand_charge["IAS_C4"] == -1
    # Sphere-0 AA column stays 0 under fuzzy mode; charge is in the ligand map.
    assert spheres.compute_charge([{ias}], structure, {}, center)[0] == 0


def test_ias_free_sidechain_gets_oxt_and_carboxylate():
    structure, ias = _build_ias_free_sidechain()
    tree = NeighborSearch(list(structure.get_atoms()))
    assert not spheres.ias_cg_isopeptide_linked(ias, tree)
    # OXT -1 + OD1/OD2 carboxylate -1
    assert spheres.hetero_residue_formal_charge(ias, tree, set()) == -2


def test_protonated_oxt_not_charged():
    structure, ias = _build_ias_free_sidechain()
    # Add hydroxyl H on OXT → saturated, no -1 from OXT
    _add_atom(ias, "HXT", "H", (4.5, 0.6, 0.0), 99)
    tree = NeighborSearch(list(structure.get_atoms()))
    # Still free sidechain -1; OXT saturated → total -1
    assert spheres.hetero_residue_formal_charge(ias, tree, set()) == -1
