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


def _add_sugar_phosphate(res, origin, serial, with_p=True):
    """Minimal sugar (+ optional phosphate) at ``origin`` (P position)."""
    o = np.asarray(origin, dtype=float)
    if with_p:
        serial = _add_atom(res, "P", "P", o + (0.0, 0.0, 0.0), serial)
        serial = _add_atom(res, "OP1", "O", o + (1.55, 0.0, 0.0), serial)
        serial = _add_atom(res, "OP2", "O", o + (-0.78, 1.34, 0.0), serial)
        serial = _add_atom(res, "O5'", "O", o + (-0.78, -1.34, 0.0), serial)
    serial = _add_atom(res, "C5'", "C", o + (-2.2, -1.8, 0.0), serial)
    serial = _add_atom(res, "C4'", "C", o + (-3.5, -1.0, 0.0), serial)
    serial = _add_atom(res, "O4'", "O", o + (-4.5, -1.5, 0.8), serial)
    serial = _add_atom(res, "C1'", "C", o + (-5.2, -0.5, 0.0), serial)
    serial = _add_atom(res, "C3'", "C", o + (-3.5, 0.5, -0.5), serial)
    serial = _add_atom(res, "O3'", "O", o + (-3.5, 1.9, -0.5), serial)
    return serial


def test_sah_not_treated_as_nucleic():
    """Ribose ligands (SAH/SAM) must not get nucleic O3' caps over existing H3'."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("A")
    structure.add(model)
    model.add(chain)
    sah = Residue(("H_SAH", 801, " "), "SAH", " ")
    serial = _add_sugar_phosphate(sah, (0.0, 0.0, 0.0), 1, with_p=False)
    serial = _add_atom(sah, "H3'", "H", np.array([-3.5, 2.7, -0.5]), serial)
    chain.add(sah)

    assert not spheres.is_polymer_nucleotide(sah)
    spheres.cap_chains(model, {sah}, capping=1, ligand_charge={})
    assert not spheres.has_res_atom(sah, "HO3'")


def test_a1c_is_polymer_nucleotide_with_phosphate_charge():
    """Modified base A1C (sugar + P) is nucleic; phosphate -1; no protein N-term."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("B")
    structure.add(model)
    model.add(chain)
    a1c = Residue((" ", 48, " "), "A1C", " ")
    serial = _add_sugar_phosphate(a1c, (0.0, 0.0, 0.0), 1, with_p=True)
    # Base N so a naive protein N-terminus rule would fire if not gated
    serial = _add_atom(a1c, "N1", "N", np.array([-6.5, -0.5, 0.0]), serial)
    chain.add(a1c)

    assert spheres.is_polymer_nucleotide(a1c)
    center = spheres.CenterResidue("exact:A1C_B48")
    charge = spheres.compute_charge([{a1c}], structure, {}, center)
    assert charge == [-1]


def test_strand_u_a1c_g_caps():
    """U–A1C–G fragment: HP only on 5' U; HO3' only on 3' G; internal uncapped."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("B")
    structure.add(model)
    model.add(chain)

    residues = []
    specs = [
        ("U", 47, (0.0, 0.0, 0.0)),
        ("A1C", 48, (10.0, 0.0, 0.0)),
        ("G", 49, (20.0, 0.0, 0.0)),
    ]
    serial = 1
    for name, seq, origin in specs:
        res = Residue((" ", seq, " "), name, " ")
        serial = _add_sugar_phosphate(res, origin, serial, with_p=True)
        chain.add(res)
        residues.append(res)

    # Link phosphodiesters: move downstream P next to upstream O3' along +y
    # so OP1/OP2 stay outside the upstream P's 2.0 Å valence sphere.
    for up, down in zip(residues, residues[1:]):
        o3 = spheres.get_res_atom(up, "O3'").get_coord()
        target = o3 + np.array([0.0, 1.6, 0.0])
        delta = target - down["P"].get_coord()
        for atom_name in ("P", "OP1", "OP2", "O5'"):
            down[atom_name].set_coord(down[atom_name].get_coord() + delta)

    cluster = set(residues)
    spheres.cap_chains(model, cluster, capping=1, ligand_charge={})
    u47, a1c, g49 = residues
    assert u47.has_id("HP"), "5' phosphate should receive P–H cap"
    assert spheres.has_res_atom(g49, "HO3'"), "3' O3' should receive H cap"
    assert not a1c.has_id("HP")
    assert not spheres.has_res_atom(a1c, "HO3'")
    assert not spheres.has_res_atom(u47, "HO3'")


def test_cys_a1c_thioether_skips_thiolate_charge():
    """CYS SG–A1C C6 adduct is neutral thioether (no RGP needed)."""
    structure = Structure("test")
    model = Model(0)
    chain_a = Chain("A")
    chain_b = Chain("B")
    structure.add(model)
    model.add(chain_a)
    model.add(chain_b)

    cys = Residue((" ", 321, " "), "CYS", " ")
    serial = 1
    serial = _add_atom(cys, "N", "N", (0.0, 0.0, 0.0), serial)
    # Amide H so missing-backbone-H does not add an extra −1 (see RGP tests).
    serial = _add_atom(cys, "H", "H", (0.0, 1.0, 0.0), serial)
    serial = _add_atom(cys, "CA", "C", (1.5, 0.0, 0.0), serial)
    serial = _add_atom(cys, "C", "C", (2.2, 1.2, 0.0), serial)
    serial = _add_atom(cys, "O", "O", (2.2, 2.4, 0.0), serial)
    serial = _add_atom(cys, "CB", "C", (1.5, -1.5, 0.0), serial)
    serial = _add_atom(cys, "SG", "S", (1.5, -3.0, 0.0), serial)
    chain_a.add(cys)

    a1c = Residue((" ", 48, " "), "A1C", " ")
    serial = _add_sugar_phosphate(a1c, (10.0, 0.0, 0.0), serial, with_p=True)
    # C6 within covalent distance of SG (~1.84 A)
    serial = _add_atom(a1c, "C6", "C", (1.5, -4.84, 0.0), serial)
    chain_b.add(a1c)

    center = spheres.CenterResidue("exact:CYS_A321-A1C_B48")
    charge = spheres.compute_charge(
        [{cys, a1c}], structure, {}, center, residues={cys, a1c}
    )
    # A1C phosphate -1; CYS thioether 0 (not thiolate -1)
    assert charge == [-1]
    from Bio.PDB.NeighborSearch import NeighborSearch

    assert spheres.check_cys_nucleotide_thioether(
        cys, NeighborSearch(list(structure.get_atoms()))
    )


def test_free_cys_still_thiolate_without_nucleotide():
    structure = Structure("test")
    model = Model(0)
    chain = Chain("A")
    structure.add(model)
    model.add(chain)
    cys = Residue((" ", 271, " "), "CYS", " ")
    serial = 1
    serial = _add_atom(cys, "N", "N", (0.0, 0.0, 0.0), serial)
    serial = _add_atom(cys, "H", "H", (0.0, 1.0, 0.0), serial)
    serial = _add_atom(cys, "CA", "C", (1.5, 0.0, 0.0), serial)
    serial = _add_atom(cys, "C", "C", (2.2, 1.2, 0.0), serial)
    serial = _add_atom(cys, "O", "O", (2.2, 2.4, 0.0), serial)
    serial = _add_atom(cys, "CB", "C", (1.5, -1.5, 0.0), serial)
    serial = _add_atom(cys, "SG", "S", (1.5, -3.0, 0.0), serial)
    chain.add(cys)

    center = spheres.CenterResidue("exact:CYS_A271")
    charge = spheres.compute_charge([{cys}], structure, {}, center)
    assert charge == [-1]
