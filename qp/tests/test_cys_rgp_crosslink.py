"""CYS–ligand covalent crosslink charge via Protoss RGP (R#) matching."""

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


def _build_cys_ligand_pair():
    """CYS without HG; SG covalently near ligand carbon (no disulfide)."""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("A")
    structure.add(model)
    model.add(chain)

    cys = Residue((" ", 186, " "), "CYS", " ")
    serial = 1
    serial = _add_atom(cys, "N", "N", (0.0, 0.0, 0.0), serial)
    serial = _add_atom(cys, "H", "H", (0.0, 1.0, 0.0), serial)
    serial = _add_atom(cys, "CA", "C", (1.5, 0.0, 0.0), serial)
    serial = _add_atom(cys, "C", "C", (2.2, 1.2, 0.0), serial)
    serial = _add_atom(cys, "O", "O", (2.2, 2.4, 0.0), serial)
    serial = _add_atom(cys, "CB", "C", (1.5, -1.5, 0.0), serial)
    serial = _add_atom(cys, "SG", "S", (1.5, -3.0, 0.0), serial)
    chain.add(cys)

    lig = Residue(("H_8NR", 505, " "), "8NR", " ")
    serial = _add_atom(lig, "CAJ", "C", (1.5, -4.7, 0.0), serial)
    serial = _add_atom(lig, "C1", "C", (1.5, -6.2, 0.0), serial)
    chain.add(lig)
    return structure, cys, lig


def test_free_cys_without_hg_is_thiolate_minus_one():
    structure, cys, lig = _build_cys_ligand_pair()
    center = spheres.CenterResidue("exact:CYS_A186")
    # Ligand present but no RGP metadata → CYS still counted as thiolate.
    charge = spheres.compute_charge(
        [{cys, lig}], structure, {"8NR_A505": 0}, center, RGP_atoms={}
    )
    assert charge == [-1]


def test_rgp_linked_cys_skips_thiolate_charge():
    structure, cys, lig = _build_cys_ligand_pair()
    center = spheres.CenterResidue("exact:CYS_A186-8NR_A505")
    rgp = {
        "8NR_A505": {
            50: {
                "atom": cys["SG"],
                "linking_atom": lig["CAJ"],
                "coord": cys["SG"].get_coord(),
                "linking_atom_coord": lig["CAJ"].get_coord(),
            }
        }
    }
    charge = spheres.compute_charge(
        [{cys, lig}],
        structure,
        {"8NR_A505": 0},
        center,
        residues={cys, lig},
        RGP_atoms=rgp,
    )
    assert charge == [0]


def test_rgp_ignored_when_ligand_absent_from_cluster():
    structure, cys, lig = _build_cys_ligand_pair()
    center = spheres.CenterResidue("exact:CYS_A186")
    rgp = {
        "8NR_A505": {
            50: {
                "atom": cys["SG"],
                "linking_atom": lig["CAJ"],
            }
        }
    }
    # Ligand not in residues/spheres → RGP_flag requires name in res_keys.
    charge = spheres.compute_charge(
        [{cys}],
        structure,
        {"8NR_A505": 0},
        center,
        residues={cys},
        RGP_atoms=rgp,
    )
    assert charge == [-1]
