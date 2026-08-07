import pytest
import numpy as np
import os
import glob
import filecmp

from qp.cluster import spheres


def check_clusters(path, out, metal_ids):
    for metal in metal_ids:
        assert os.path.isdir(os.path.join(path, metal)), f"Cluster center {metal} not found"
        sphere_count = len(glob.glob(os.path.join(path, metal, "?.pdb")))
        for i in range(sphere_count):
            expected_pdb = os.path.join(path, metal, f"{i}.pdb")
            output_pdb = os.path.join(out, metal, f"{i}.pdb")
            assert filecmp.cmp(expected_pdb, output_pdb), f"Sphere {i} PDB does not match expected"

    # Expected charge.csv includes additional ligands, which must be excluded
    expected_charge = os.path.join(path, "charge.csv")
    output_charge = os.path.join(out, "charge.csv")
    with open(expected_charge, "r") as e, open(output_charge, "r") as o:
        expected_lines = e.readlines()
        output_lines = o.readlines()
        assert expected_lines[:len(output_lines)] == output_lines, "Charge does not match expected"

    expected_count = os.path.join(path, "count.csv")
    output_count = os.path.join(out, "count.csv")
    assert filecmp.cmp(expected_count, output_count), "Residue count does not match expected"


def residue_atom_lines(pdb_path, resname, chain, resnum):
    """Return the ATOM/HETATM lines for one residue in a PDB file."""
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if (
                line.startswith(("ATOM", "HETATM"))
                and line[17:20].strip() == resname
                and line[21].strip() == chain
                and int(line[22:26]) == resnum
            ):
                lines.append(line)
    return lines


@pytest.mark.parametrize("sample_cluster", [
    ("1sp9", ("A446", "B446")),
    ("2q4a", ("A901", "B902")),
    ("3a8g", ("A301",))
], indirect=True)
def test_extract_clusters(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"],
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("4ilv", ("A301", "B301"))], indirect=True)
def test_cap_heavy(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], capping=2, 
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("1lm6", ("A204",))], indirect=True)
def test_box_plot(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"],
        smooth_method="box_plot"
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2r6s", ("A501",))], indirect=True)
def test_dbscan(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 
        smooth_method="dbscan", eps=6, min_samples=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [
    ("2chb", (
        "A1_A2_A3_A4",
        "B1_B2_B3_B4_B5",
        "C1_C2_C3_C4",
        "I1_I2_I3_I4_I5",
        "J1_J2_J3_J4_J5"
    )),
    ("4z42", ("C601_C602", "F601_F602", "I601_I602", "L601_L602"))
], indirect=True)
def test_merge_centers(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["BGC", "GAL", "NGA", "SIA", "NI"],
        merge_cutoff=4.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2fd8", ("B501_B502_B503",))], indirect=True)
def test_prune_atoms(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["DT", "MA7"],
        max_atom_count=102, merge_cutoff=2.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


def test_force_include_residues(tmpdir):
    """force_include_residues should force-include a residue outside the
    default spheres, cap it like any other extracted residue, and protect
    it from max_atom_count pruning.

    Note: unlike the other tests in this file, this constructs a real
    CenterResidue (the current extract_clusters API) rather than passing a
    plain resname list.
    """
    from qp.cluster.spheres import CenterResidue

    path = os.path.join(os.path.dirname(__file__), "samples", "3a8g")
    pdb_path = os.path.join(path, "Protoss", "3a8g_protoss.pdb")
    center_residue = CenterResidue("FE")

    def cluster_pdb(out):
        return os.path.join(out, "A301", "A301.pdb")

    # GLU_A9 is not reached by sphere growth (confirmed against the golden
    # 3a8g/A301 spheres, which only cover residues near the active site)
    baseline = tmpdir.mkdir("baseline")
    spheres.extract_clusters(
        pdb_path, str(baseline), center_residue,
        smooth_method="dummy_atom", mean_distance=3
    )
    assert not residue_atom_lines(cluster_pdb(str(baseline)), "GLU", "A", 9)

    # force_include_residues should force it in and cap it
    added = tmpdir.mkdir("added")
    spheres.extract_clusters(
        pdb_path, str(added), center_residue,
        smooth_method="dummy_atom", mean_distance=3,
        force_include_residues=["GLU_A9"]
    )
    added_lines = residue_atom_lines(cluster_pdb(str(added)), "GLU", "A", 9)
    source_lines = residue_atom_lines(pdb_path, "GLU", "A", 9)
    assert added_lines, "GLU_A9 was not added to the cluster"
    assert len(added_lines) > len(source_lines), "GLU_A9 was not capped"

    # Protected from max_atom_count pruning, even though it's one of the
    # most distant residues from the active site and would normally be
    # pruned first
    pruned = tmpdir.mkdir("pruned")
    spheres.extract_clusters(
        pdb_path, str(pruned), center_residue,
        smooth_method="dummy_atom", mean_distance=3,
        force_include_residues=["GLU_A9"], max_atom_count=60
    )
    assert residue_atom_lines(cluster_pdb(str(pruned)), "GLU", "A", 9), (
        "GLU_A9 was pruned despite being explicitly requested"
    )


def test_force_include_residues_multiple(tmpdir):
    """force_include_residues should support force-including more than one
    residue at once, spanning different chains."""
    from qp.cluster.spheres import CenterResidue

    path = os.path.join(os.path.dirname(__file__), "samples", "3a8g")
    pdb_path = os.path.join(path, "Protoss", "3a8g_protoss.pdb")
    center_residue = CenterResidue("FE")

    def cluster_pdb(out):
        return os.path.join(out, "A301", "A301.pdb")

    # Both are absent from the default spheres (same reasoning as
    # test_force_include_residues), one on each chain
    residues_to_add = [("GLU", "A", 9), ("GLY", "B", 3)]

    added = tmpdir.mkdir("added")
    spheres.extract_clusters(
        pdb_path, str(added), center_residue,
        smooth_method="dummy_atom", mean_distance=3,
        force_include_residues=[f"{resname}_{chain}{resnum}" for resname, chain, resnum in residues_to_add]
    )

    for resname, chain, resnum in residues_to_add:
        added_lines = residue_atom_lines(cluster_pdb(str(added)), resname, chain, resnum)
        source_lines = residue_atom_lines(pdb_path, resname, chain, resnum)
        assert added_lines, f"{resname}_{chain}{resnum} was not added to the cluster"
        assert len(added_lines) > len(source_lines), f"{resname}_{chain}{resnum} was not capped"


def test_force_include_residues_already_present(tmpdir):
    """force_include_residues should not write a residue twice when it's
    already part of the cluster via normal sphere growth."""
    from qp.cluster.spheres import CenterResidue

    path = os.path.join(os.path.dirname(__file__), "samples", "3a8g")
    pdb_path = os.path.join(path, "Protoss", "3a8g_protoss.pdb")
    center_residue = CenterResidue("FE")

    def cluster_pdb(out):
        return os.path.join(out, "A301", "A301.pdb")

    # ALA_A113 is already part of the default sphere 1 (confirmed against
    # the golden 3a8g/A301/1.pdb), flanked on both sides by residues that
    # are also already in the cluster, so it isn't even capped -- any
    # atom-count change here can only come from duplicate writing
    baseline = tmpdir.mkdir("baseline")
    spheres.extract_clusters(
        pdb_path, str(baseline), center_residue,
        smooth_method="dummy_atom", mean_distance=3
    )
    baseline_lines = residue_atom_lines(cluster_pdb(str(baseline)), "ALA", "A", 113)
    assert baseline_lines, "test setup assumption broken: ALA_A113 should already be in the default cluster"

    forced = tmpdir.mkdir("forced")
    spheres.extract_clusters(
        pdb_path, str(forced), center_residue,
        smooth_method="dummy_atom", mean_distance=3,
        force_include_residues=["ALA_A113"]
    )
    forced_lines = residue_atom_lines(cluster_pdb(str(forced)), "ALA", "A", 113)

    assert len(forced_lines) == len(baseline_lines), (
        "ALA_A113 was written a different number of times when force-included "
        "despite already being part of the cluster (expected no duplication)"
    )


@pytest.mark.parametrize("sample_cluster", [
    ("2chb", (
        "A1_A2_A3_A4",
        "B1_B2_B3_B4_B5",
        "C1_C2_C3_C4",
        "I1_I2_I3_I4_I5",
        "J1_J2_J3_J4_J5"
    ))
], indirect=True)
def test_cluster_name_template(tmpdir, sample_cluster):
    """cluster_name_template should produce short, collision-free names
    (avoiding the long residue-concatenated names in the default case,
    which can exceed path length limits when many centers are merged),
    while the underlying sphere geometry stays identical."""
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, ["BGC", "GAL", "NGA", "SIA", "NI"],
        merge_cutoff=4.0,
        smooth_method="dummy_atom", mean_distance=3,
        first_sphere_radius=4.0,
        cluster_name_template="A_{radius}"
    )

    expected_names = ["A_4", "A_4_1", "A_4_2", "A_4_3", "A_4_4"]
    for name in expected_names:
        assert os.path.isdir(os.path.join(tmpdir, name)), f"Expected cluster dir {name} not found"

    # No unexpected directories, and no collisions/overwrites occurred
    generated_dirs = sorted(
        d for d in os.listdir(tmpdir)
        if os.path.isdir(os.path.join(tmpdir, d))
    )
    assert generated_dirs == sorted(expected_names)

    # Sphere geometry should be identical to the default-naming output, just
    # under new directory names. Centers aren't guaranteed to be processed
    # in a fixed order (get_center_residues iterates over a set), so match
    # each generated cluster to its expected counterpart by content rather
    # than by position.
    remaining_expected = list(metal_ids)
    for new_name in expected_names:
        output_pdb_0 = os.path.join(tmpdir, new_name, "0.pdb")
        match = None
        for old_name in remaining_expected:
            expected_pdb_0 = os.path.join(path, old_name, "0.pdb")
            if filecmp.cmp(expected_pdb_0, output_pdb_0):
                match = old_name
                break
        assert match is not None, f"Generated cluster {new_name} does not match any expected cluster"
        remaining_expected.remove(match)

        sphere_count = len(glob.glob(os.path.join(path, match, "?.pdb")))
        for i in range(sphere_count):
            expected_pdb = os.path.join(path, match, f"{i}.pdb")
            output_pdb = os.path.join(tmpdir, new_name, f"{i}.pdb")
            assert filecmp.cmp(expected_pdb, output_pdb), (
                f"Sphere {i} PDB for {match} -> {new_name} does not match expected"
            )
    assert not remaining_expected, "Not all expected clusters were matched"


def test_cluster_name_template_bad_field(tmpdir):
    """An unknown field in the template should raise a clear error rather
    than failing silently or with a confusing traceback."""
    with pytest.raises(ValueError, match="cluster_name_template"):
        spheres.extract_clusters(
            os.path.join(os.path.dirname(__file__), "samples", "1lm6", "Protoss", "1lm6_protoss.pdb"),
            tmpdir, ["FE", "FE2"],
            smooth_method="box_plot",
            cluster_name_template="A_{not_a_real_field}"
        )
