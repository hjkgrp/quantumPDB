"""Tests for temporary mmCIF → PDB conversion."""

import json
import os
from unittest.mock import MagicMock

import pytest
from Bio.PDB import PDBParser

from qp.structure import setup
from qp.structure.mmcif_to_pdb import (
    OversizedStructureError,
    convert_mmcif_to_pdb,
    remap_sidecar_path,
)

# Example lives at the NeuralBioChem workspace root (sibling of src/)
CIF_21ZQ = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "21ZQ.cif")
)


@pytest.mark.skipif(not os.path.isfile(CIF_21ZQ), reason="21ZQ.cif example not found")
def test_convert_21zq_mmcif(tmpdir):
    out = os.path.join(tmpdir, "21ZQ.pdb")
    info = convert_mmcif_to_pdb(CIF_21ZQ, out)

    assert os.path.isfile(out)
    assert os.path.getsize(out) > 0
    assert "A1E3R" in info["resname_map"]
    assert len(info["resname_map"]["A1E3R"]) == 3

    sidecar = remap_sidecar_path(out)
    assert os.path.isfile(sidecar)
    with open(sidecar) as handle:
        sidecar_data = json.load(handle)
    assert sidecar_data["resname_map"]["A1E3R"] == info["resname_map"]["A1E3R"]

    structure = PDBParser(QUIET=True).get_structure("21ZQ", out)
    n_atoms = sum(1 for _ in structure.get_atoms())
    assert n_atoms == 4168

    for residue in structure.get_residues():
        assert len(residue.resname) <= 3
    for model in structure:
        for chain in model:
            assert len(chain.id) == 1

    # Occupancy / B-factor columns must not be glued together
    with open(out) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                assert len(line) >= 66
                assert line[54:66].count(".") <= 2


def test_get_pdbs_accepts_cif(tmpdir):
    cif_path = os.path.join(tmpdir, "demo.cif")
    with open(cif_path, "w") as handle:
        handle.write("data_demo\n")
    output = os.path.join(tmpdir, "out")
    pdb_all = setup.get_pdbs([cif_path], output)
    assert len(pdb_all) == 1
    pdb_id, pdb_path, source_cif = pdb_all[0]
    assert pdb_id == "demo"
    assert pdb_path == os.path.join(output, "demo", "demo.pdb")
    assert source_cif == os.path.abspath(cif_path)


def test_ensure_structure_pdb_converts_local_cif(tmpdir):
    if not os.path.isfile(CIF_21ZQ):
        pytest.skip("21ZQ.cif example not found")

    pdb_path = os.path.join(tmpdir, "21ZQ", "21ZQ.pdb")
    status = setup.ensure_structure_pdb("21ZQ", pdb_path, source_cif=CIF_21ZQ)
    assert status == "converted"
    assert os.path.isfile(pdb_path)
    assert setup.ensure_structure_pdb("21ZQ", pdb_path, source_cif=CIF_21ZQ) == "exists"


def test_fetch_pdb_mmcif_fallback(tmpdir, monkeypatch):
    """When classic PDB 404s, fetch mmCIF and convert."""
    if not os.path.isfile(CIF_21ZQ):
        pytest.skip("21ZQ.cif example not found")

    with open(CIF_21ZQ) as handle:
        cif_text = handle.read()

    def fake_get(url, timeout=30):
        response = MagicMock()
        if url.endswith(".pdb"):
            response.status_code = 404
            response.text = ""
        elif url.endswith(".cif"):
            response.status_code = 200
            response.text = cif_text
        else:
            response.status_code = 500
            response.text = ""
        return response

    monkeypatch.setattr(setup.requests, "get", fake_get)
    out = os.path.join(tmpdir, "21ZQ", "21ZQ.pdb")
    status = setup.fetch_pdb("21ZQ", out)
    assert status == "mmcif"
    assert os.path.isfile(out)
    assert os.path.isfile(os.path.join(tmpdir, "21ZQ", "21ZQ.cif"))
    structure = PDBParser(QUIET=True).get_structure("21ZQ", out)
    assert sum(1 for _ in structure.get_atoms()) == 4168


def test_fetch_pdb_invalid_id(tmpdir, monkeypatch):
    def fake_get(url, timeout=30):
        response = MagicMock()
        response.status_code = 404
        response.text = ""
        return response

    monkeypatch.setattr(setup.requests, "get", fake_get)
    out = os.path.join(tmpdir, "XXXX.pdb")
    with pytest.raises(ValueError, match="XXXX"):
        setup.fetch_pdb("XXXX", out)


def test_oversized_structure_raises(tmpdir, monkeypatch):
    """Oversized conversions raise OversizedStructureError for batch skip handling."""
    if not os.path.isfile(CIF_21ZQ):
        pytest.skip("21ZQ.cif example not found")

    monkeypatch.setattr("qp.structure.mmcif_to_pdb.MAX_PDB_ATOMS", 1)
    out = os.path.join(tmpdir, "too_big.pdb")
    with pytest.raises(OversizedStructureError, match="atoms"):
        convert_mmcif_to_pdb(CIF_21ZQ, out)
    assert not os.path.isfile(out)
