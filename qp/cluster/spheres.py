"""Extract active site coordination sphere clusters

**Usage**::

    >>> from qp.cluster import spheres

    >>> spheres.extract_clusters(
    ...     "path/to/PDB.pdb", 
    ...     "path/to/out/dir/", 
    ...     center_residues=["FE", "FE2"], # List of resnames of the residues to use as the cluster center
    ...     sphere_count=2,              # Number of spheres to extract
    ...     ligands=["AKG"]       # PDB IDs of additional ligands
    ... )

Extracting clusters leaves open valences in the outermost sphere. Capping may be
performed by specifying ``capping`` in ``spheres.extract_clusters``:

* 0. No capping. (Default)
* 1. Cap with hydrogens. Proteins: N/C termini. Non-ligand A/C/G/U: O3'–H
  at strand breaks.
* 2. Cap with ACE/NME groups (proteins only). Non-ligand A/C/G/U still use
  O3'–H caps when ``capping`` is non-zero.
"""

import os
from functools import reduce
from typing import Set, Literal, Optional, List, Dict, Any, Tuple
import numpy as np
from Bio.PDB import PDBParser, Polypeptide, PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.NeighborSearch import NeighborSearch
from scipy.spatial import Voronoi
from qp.cluster import struct_to_file
from sklearn.cluster import DBSCAN
from qp.structure.mmcif_to_pdb import expand_resnames_for_matching, normalize_center_key


RANDOM_SEED = 66265
HX_BOND_LENGTH = {
    "C": 1.09,
    "N": 1.00,
    "O": 0.98,
    "S": 1.35,
}

# RNA bases that may appear as polymer residues without a Protoss ligand entry.
RNA_POLYMER_RESNAMES = frozenset({"A", "C", "G", "U"})
# m7GTP / similar residues: triphosphate + N7-methyl formal charges.
MGT_RESNAMES = frozenset({"MGT"})

O_H_BOND = 0.97
_TET_COS = np.cos(np.deg2rad(109.5))
_TET_SIN = np.sin(np.deg2rad(109.5))

CHARGE_DEBUG_FLAG = False


def charge_debug(msg, res=None):
    """Print charge-assignment debug messages when ``CHARGE_DEBUG_FLAG`` is set."""
    if CHARGE_DEBUG_FLAG:
        print(msg, "" if res is None else make_res_key(res))

class CenterResidue:
    def __init__(self, center_residue: str, resname_map=None):
        """Parse a center residue definition string.

        Matching mode is chosen as follows:

        * ``exact:FE_A199`` --- **strict** single-residue selection. The
          ``exact:`` prefix is required for a single ``RESNAME_CHAINID`` key
          because bare ``FE_A199`` would otherwise be ambiguous with fuzzy
          resnames (e.g. CCD code ``A199``).
        * ``CU_A357-CU_A358`` --- **strict** multi-residue list. Dash-separated
          exact keys do not need the ``exact:`` prefix.
        * ``exact:FE_A155-HIS_A93`` --- also strict; the prefix is allowed on
          multi-residue lists too.
        * ``FE`` / ``FE_FE2`` --- **fuzzy** mode: underscore-separated residue
          names matched against HETATM records only.

        When ``resname_map`` is provided (from an mmCIF→PDB remap sidecar),
        original longer residue names (e.g. 5-letter CCD codes) are also
        accepted and matched against the remapped 3-letter names present in
        the structure. Other pipeline stages that use resnames are unchanged.

        Parameters
        ----------
        center_residue : str
            Center definition string from the config or CSV input.
        resname_map : dict, optional
            Mapping of original residue names → PDB residue names
            (``{original: mapped}``), typically from ``*_mmcif_remap.json``.
        """
        self.center_residue_str = center_residue
        self.resname_map = dict(resname_map or {})

        raw = center_residue.strip()
        exact_prefix = raw.lower().startswith("exact:")
        if exact_prefix:
            raw = raw.split(":", 1)[1].strip()
            if not raw:
                raise ValueError(
                    "Center residue 'exact:' prefix must be followed by a "
                    "RESNAME_CHAINID key (e.g. 'exact:FE_A199')."
                )

        residue_list = raw.split("-")
        if exact_prefix or len(residue_list) > 1:
            self.mode = "strict"
            self.residue_list = residue_list
        else:
            self.mode = "fuzzy"
            self.residue_list = raw.split("_")

    def __str__(self):
        return self.center_residue_str
    
    def __repr__(self):
        return self.center_residue_str

    def _structure_resnames(self):
        """Residue names that should match structure ``get_resname()`` values."""
        return expand_resnames_for_matching(self.residue_list, self.resname_map)

    def _normalize_strict_key(self, key: str) -> str:
        """Rewrite a strict center key through ``resname_map`` when needed."""
        return normalize_center_key(key, self.resname_map)

    def __contains__(self, res: Residue):
        """Test whether a residue matches this center definition.

        In fuzzy mode, matches any HETATM residue whose name appears in the
        residue list (or maps to that name via ``resname_map``). In strict
        mode, matches by exact ``RESNAME_CHAINID`` key, accepting original
        mmCIF residue names when a remap is available.

        Parameters
        ----------
        res : Bio.PDB.Residue.Residue
            Residue to test.

        Returns
        -------
        bool
        """
        if self.mode == "fuzzy":
            return res.get_resname() in self._structure_resnames() and res.id[0] != ' '
        elif self.mode == "strict":
            key = make_res_key(res)
            allowed = {self._normalize_strict_key(tok) for tok in self.residue_list}
            return key in allowed


def get_grid_coord_idx(coord, coord_min, mean_distance):
    """
    Compute a point's position in a 1D grid
    
    Parameters
    ----------
    coord: float
        One coordinate of a point
    coord_min: float
        The minimum coordinate of the grid
    mean_distance: float
        The distance between neighbors in the grid

    Returns
    -------
    idx: int
        The idx of the point in the grid
    """
    return int((coord - coord_min + mean_distance * 0.5) // mean_distance)


def get_grid(coords, mean_distance):
    """
    Compute the grid's parameter for a given 1D point list

    Parameters
    ----------
    coords: numpy.array
        The coordinates of the 1D point list
    mean_distance: float
        The distance between neighbors in the grid
    
    Returns
    -------
    coord_min: float
        The minimum coordinate of the grid
    coord_max: float
        The maximum coordinate of the grid
    grid: numpy.array
        The 1D grid
    """
    coord_min, coord_max = coords.min() - mean_distance, coords.max() + mean_distance
    npoints = get_grid_coord_idx(coord_max, coord_min, mean_distance)
    return coord_min, coord_max, coord_min + np.linspace(0, npoints - 1, npoints) * mean_distance


def visualize_dummy(dummy, path):
    """Write dummy atom positions to an XYZ file for debugging.

    Parameters
    ----------
    dummy : list of array-like
        3D coordinates of dummy atoms.
    path : str
        Directory where ``dummy.xyz`` will be written.
    """
    with open(path + "/dummy.xyz", "w") as f:
        f.write(f"{len(dummy)}\n\n")
        for coord in dummy:
            f.write(f"He {coord[0]} {coord[1]} {coord[2]}\n")


def fill_dummy(points, mean_distance=3, noise_amp=0.2):
    """
    Fill dummy atoms in a point cloud

    Parameters
    ----------
    points: numpy.array
        The 3D coordinates of the point cloud
    mean_distance: float
        The distance between neighbors in the grid
    noise_amp: float
        The amplitude of the noise of dummy atoms' position
    
    Returns
    -------
    points:
        The 3D coordinates of the point cloud filled with dummy atoms
    """
    conf = np.stack(points, axis=0)
    x_min, _, x_grids = get_grid(conf[:,0], mean_distance)
    y_min, _, y_grids = get_grid(conf[:,1], mean_distance)
    z_min, _, z_grids = get_grid(conf[:,2], mean_distance)
    flags = np.ones((len(x_grids), len(y_grids), len(z_grids)), dtype=bool)
    for point in points:
        x, y, z = point
        flags[
            get_grid_coord_idx(x, x_min, mean_distance),
            get_grid_coord_idx(y, y_min, mean_distance), 
            get_grid_coord_idx(z, z_min, mean_distance)
        ] = False
    flags = flags.flatten()
    np.random.seed(RANDOM_SEED)
    noise = mean_distance * noise_amp * (np.random.rand(len(x_grids), len(y_grids), len(z_grids), 3) - 0.5)
    dummy = np.stack(np.meshgrid(x_grids, y_grids, z_grids, indexing="ij"), axis=-1)
    dummy = (dummy + noise).reshape(-1, 3)
    dummy = dummy[flags, :]
    return np.concatenate([points, dummy], axis=0)


def calc_dist(point_a, point_b):
    """
    Calculate the Euclidean distance between two points

    Parameters
    ----------
    point_a: numpy.array
        Point A
    point_b: numpy.array
        Point B

    Returns
    -------
    dist: float
        The Euclidean distance between Point A and Point B    

    """
    return np.linalg.norm(point_a - point_b)


def voronoi(model, center_residue: CenterResidue, ligands, smooth_method, output_path, **smooth_params):
    """Compute the Voronoi tessellation of a protein structure.

    Builds a Voronoi diagram from all atoms in the model and returns an
    adjacency list of neighboring atoms. Optionally applies smoothing via
    dummy atoms to improve boundary handling for surface-exposed residues.

    Parameters
    ----------
    model : Bio.PDB.Model.Model
        Protein structure model containing all residues.
    center_residue : CenterResidue
        Center residue definition (currently unused but kept for API).
    ligands : list
        Ligand residue names to include (currently unused but kept for API).
    smooth_method : str
        Smoothing method: ``'dummy_atom'`` fills voids with dummy atoms,
        otherwise no smoothing is applied.
    output_path : str
        Directory for debug output (e.g., ``dummy.xyz``).
    **smooth_params
        Additional parameters passed to :func:`fill_dummy` (e.g.,
        ``mean_distance``, ``noise_amp``).

    Returns
    -------
    dict
        Adjacency list mapping each ``Bio.PDB.Atom`` to a list of
        ``(neighbor_atom, distance)`` tuples for all Voronoi neighbors.
    """
    atoms = []
    points = []
    for res in model.get_residues():
        # if Polypeptide.is_aa(res) or \
        #     res.get_resname() in center_residues or \
        #     (res.get_resname() in ligands):
        for atom in res.get_unpacked_list(): # includes atoms from multiple conformations
            atoms.append(atom)
            points.append(atom.get_coord())

    points_count = len(points)
    if smooth_method == "dummy_atom":
        new_points = fill_dummy(points, **smooth_params)
        visualize_dummy(new_points, output_path)
        vor = Voronoi(new_points)
        # vor = Voronoi(points)
    else:
        vor = Voronoi(points)

    # Plot Voronoi diagrams (2D and 3D)
    # plot_voronoi_2d(new_points, points_count, output_path)  # Plot 2D Voronoi excluding dummy atoms

    neighbors = {}
    for a, b in vor.ridge_points:
        if a < points_count and b < points_count:
            dist = calc_dist(points[a], points[b])
            neighbors.setdefault(atoms[a], []).append((atoms[b], dist))
            neighbors.setdefault(atoms[b], []).append((atoms[a], dist))
    return neighbors


def merge_centers(cur, search, seen, radius=0.0):
    """Recursively merge center residues that are within a distance cutoff.

    Starting from ``cur``, finds all neighboring residues within ``radius``
    using a NeighborSearch and recursively merges them into a single center
    set. Hydrogen atoms are excluded from the distance search.

    Parameters
    ----------
    cur : Bio.PDB.Residue.Residue
        Starting residue.
    search : Bio.PDB.NeighborSearch
        Spatial search object built from center atom coordinates.
    seen : set
        Already-visited residues (modified in place to prevent cycles).
    radius : float, optional
        Merge distance cutoff in angstroms (default 0.0 disables merging).

    Returns
    -------
    set
        Set of residues that form the merged center.
    """
    if radius == 0.0:
        return {cur}

    center = {cur}
    seen.add(cur)
    nxt = set()
    for atom in cur.get_unpacked_list():
        if atom.element == "H":
            continue
        for res in search.search(atom.get_coord(), radius, "R"):
            nxt.add(res)
            
    for res in nxt:
        if res not in seen:
            center |= merge_centers(res, search, seen, radius)
    return center


def get_center_residues(model, center_residue: CenterResidue, merge_cutoff=0.0):
    """Find all center residues in a model and optionally merge nearby ones.

    Scans all residues in the model for matches against the center definition,
    then groups them using :func:`merge_centers` if ``merge_cutoff > 0``.

    Parameters
    ----------
    model : Bio.PDB.Model.Model
        Protein structure model to search.
    center_residue : CenterResidue
        Center residue definition.
    merge_cutoff : float, optional
        Distance cutoff in angstroms for merging nearby centers (default 0.0).

    Returns
    -------
    list of set
        Each set contains the residues forming one (possibly merged) center.
        Returns an empty list if no matching center is found.
    """
    found = set()
    center_atoms = []
    for res in model.get_residues():
        # We will assume the ligand is labeled as a heteroatom with res.id[0] != ' '
        if res in center_residue:
            found.add(res)
            center_atoms.extend([atom for atom in res.get_unpacked_list() if atom.element != "H"])
    if not len(center_atoms):
        print("> WARNING: No matching cluster center found. Skipping cluster.")
        return []
    
    search = NeighborSearch(center_atoms)
    seen = set()
    centers = []
    for res in found:
        if res not in seen:
            centers.append(merge_centers(res, search, seen, merge_cutoff))
    return centers


def box_outlier_thres(data, coeff=1.5):
    """
    Compute the threshold for the boxplot outlier detection method

    Parameters
    ----------
    data: list
        The data for the boxplot statistics
    coeff: float
        The coefficient for the outlier criterion from quartiles of the data

    Returns
    -------
    lb: float
        the lower bound of non-outliers
    ub: float
        the upper bound of non-outliers

    """
    Q3 = np.quantile(data, 0.75)
    Q1 = np.quantile(data, 0.25)
    IQR = Q3 - Q1
    return Q1 - coeff * IQR, Q3 + coeff * IQR


def check_NC(atom, metal):
    """
    Check if a nitrogen / carbon atom in the first sphere is coordinated.

    If the nitrogen / carbon atom is in the backbone, or it's the nearest atom to the metal
    among all atoms in its residue, it's considered coordinated.

    Parameters
    ----------
    atom: Bio.PDB.Atom
        The nitrogen / carbon atom to be checked
    metal:
        The metal atom (coordination center)

    Returns
    -------
    flag: bool
        Whether the atom is coordinated or not
    """
    if atom.get_name() == "N":
        # TODO: check special coordinated nitrogens in the backbone
        return True
    else:
        ref_dist = calc_dist(atom.get_coord(), metal.get_coord())
        res = atom.get_parent()
        for atom in res.get_unpacked_list():
            if calc_dist(atom.get_coord(), metal.get_coord()) < ref_dist:
                return False
        return True


def get_next_neighbors(
    start, neighbors, sphere_count, ligands,
    first_sphere_radius=4,
    smooth_method="boxplot", 
    include_ligands=2,
    **smooth_params):
    """
    Iteratively determines spheres around a given starting atom

    Parameters
    ----------
    start: Bio.PDB.Residue
        Starting metal atom
    neighbors: dict
        Adjacency list of neighboring atoms
    sphere_count: int
        Number of spheres to extract
    ligands: list
        A list of ligands to include
    smooth_method: ("boxplot" | "dbscan" | "dummy_atom")
        The method used to smoothen the spheres
    include_ligands: int
        the mode of including ligands in the sphere
    smooth_params:
        params of the specific smooth method

    Returns
    -------
    metal_id: str
        Active site identifier
    seen: set
        Set of residues from all spheres
    spheres: list of sets
        Sets of residues separated by spheres
    """
    seen = start.copy()
    spheres = [start]
    lig_frontiers = [set()]
    lig_adds = [set()]
    start_atoms = []
    for res in start:
        start_atoms.extend(res.get_unpacked_list())
    is_metal_like = (len(start_atoms) == len(start))
    search = NeighborSearch([atom for atom in start_atoms[0].get_parent().get_parent().get_parent().get_parent().get_atoms() if atom.element != "H" and atom not in start_atoms])
    for i in range(0, sphere_count):
        # get candidate atoms in the new sphere
        nxt = set()
        lig_add = set()
        lig_frontier_atoms = set()
        if i == 0 and first_sphere_radius > 0:
            first_sphere = set()
            for center in start_atoms:
                first_sphere |= set(search.search(center=center.get_coord(), radius=first_sphere_radius, level="A"))
            for atom in first_sphere:
                if atom.get_parent() not in seen:
                    element = atom.element
                    if (
                        not is_metal_like 
                        or element in "OS" 
                        or (element in "NC" and any(check_NC(atom, center) for center in start_atoms))
                    ): # only consider coordinated atoms
                        res = atom.get_parent()
                        seen.add(res)
                        if Polypeptide.is_aa(res):
                            if include_ligands != 3 or Polypeptide.is_aa(res, standard=True):
                                nxt.add(res)
                        else:
                            if (
                                (include_ligands != 3 or res.get_resname() == "HOH") and ( # mode 3: only include center, waters, standard AAs
                                include_ligands != 1 or
                                res.get_resname() != "HOH") # mode 1: exclude all waters
                            ):
                                lig_frontier_atoms.add(atom)
                                lig_add.add(res)
        else:   
            candidates = []
            frontiers = spheres[-1] if spheres[-1] else spheres[0] # if no previous sphere, use the starting atoms
            frontier_atoms = [atom for res in frontiers for atom in res.get_unpacked_list()]
            for atom in frontier_atoms:
                for nb_atom, dist in neighbors[atom]:
                    par = nb_atom.get_parent()
                    candidates.append((nb_atom, dist, par))

            # screen candidates
            if smooth_method == "box_plot":
                dist_data = [dist for _, dist, _ in candidates]
                _, ub = box_outlier_thres(dist_data, **smooth_params)
                screened_candidates = [(atom, dist, par) for atom, dist, par in candidates if dist < ub]
            elif smooth_method == "dbscan":
                optics = DBSCAN(**smooth_params)
                X = [atom.get_coord() for atom, _, _ in candidates]
                for res in spheres[-1]:
                    for atom in res.get_unpacked_list():
                        X.append(atom.get_coord())
                cluster_idx = optics.fit_predict(X) + 1
                largest_idx = np.bincount(cluster_idx).argmax()
                screened_candidates = [(atom, dist, par) for i, (atom, dist, par) in enumerate(candidates) if cluster_idx[i] == largest_idx]
            else:
                screened_candidates = [(atom, dist, par) for atom, dist, par in candidates]

            new_ball = set()
            for atom in lig_frontiers[i]:
                new_ball |= set(search.search(center=atom.get_coord(), radius=first_sphere_radius, level="A"))
            for atom in new_ball:
                screened_candidates.append((atom, 0, atom.get_parent()))

            # build the new sphere
            new_seen = set()
            for atom, _, par in screened_candidates:
                if par not in seen:
                    new_seen.add(par)
                ## Added to include ligands in the first coordination sphere
                ## Produced unweildy clusters as some ligands are large
                ## Potentially useful in the future
                    if Polypeptide.is_aa(par):
                        nxt.add(par)
                    else:
                        if (
                            (include_ligands == 0 and par.get_resname() in ligands) or 
                            # mode 0: only include ligands in the first sphere unless specified
                            (include_ligands == 1 and par.get_resname() != "HOH") or # mode 1: exclude all waters
                            include_ligands == 2 or # mode 2: include everything
                            (include_ligands == 3 and par.get_resname() == "HOH") # mode 3: only include center, waters, standard AAs
                        ):
                            lig_frontier_atoms.add(atom)
                            lig_add.add(par)
            seen |= new_seen

        spheres.append(nxt)
        lig_adds.append(lig_add)
        lig_frontiers.append(lig_frontier_atoms)

    metal_id = []
    for res in start:
        res_id = res.get_full_id()
        chain_name = res_id[2]
        metal_index = str(res_id[3][1])
        metal_id.append(chain_name + metal_index)
    for i in range(len(spheres)):
        spheres[i] = spheres[i] | lig_adds[i]
    
    return "_".join(sorted(metal_id)), reduce(lambda x, y: x | y, spheres), spheres


def prune_atoms(center, residues, spheres, max_atom_count, ligands, kept_monomers=None):
    """Prune residues from the cluster to meet the max atom count constraint.

    Removes residues furthest from the center first, while preserving
    specified ligands and co-factors. Modifies ``residues`` and ``spheres``
    in place. Empty outer spheres are removed from the list.

    Parameters
    ----------
    center : set
        Set of central residues (used as distance reference).
    residues : set
        Set of all residues in the cluster (modified in place).
    spheres : list of set
        Residue sets by coordination sphere (modified in place).
    max_atom_count : int
        Maximum allowed total atom count in the cluster.
    ligands : list
        Ligand residue names to preserve regardless of distance.
    kept_monomers : list, optional
        Oligomer monomers that must be preserved from pruning.

    Notes
    -----
    This function operates in place and does not return a value. Residues
    are removed in order of decreasing distance from the center atoms.
    """
    if kept_monomers is None:
        kept_monomers = []

    atom_cnt = 0
    for res in residues:
        atom_cnt += len(res)
    if atom_cnt <= max_atom_count:
        return

    center_atoms = []
    for c in center:
        center_atoms.extend(c.get_unpacked_list())
    def dist(res):
        return min(atom - x for x in center_atoms for atom in res.get_unpacked_list())
                   
    prune = set()
    for res in sorted(residues, key=dist, reverse=True):
        # Check if the residue is in the ligands_to_keep list
        if res.get_resname() not in ligands and res not in kept_monomers:
            prune.add(res)
            atom_cnt -= len(res)
            if atom_cnt <= max_atom_count:
                break

    residues -= prune
    for s in spheres:
        s -= prune
    while not spheres[-1]:
        spheres.pop()


def scale_hydrogen(a, b, scale):
    """
    Replaces an atom with hydrogen, rescaling the original bond length

    Parameters
    ----------
    a: Bio.PDB.Atom
        Bonded atom to keep
    b: Bio.PDB.Atom
        Bonded atom to replace
    scale: float
        Bond length scale

    Returns
    -------
    pos: array of float
        Coordinates of new hydrogen atom
    """
    p = a.get_coord()
    q = b.get_coord()
    return scale * (q - p) + p


def get_normalized_vector(atom1: Atom, atom2: Atom) -> np.array:
    """Return the unit vector pointing from ``atom1`` to ``atom2``."""
    v = atom2.get_coord() - atom1.get_coord()
    return v / np.linalg.norm(v)


def build_hydrogen(
    parent: Residue,
    template: Optional[Residue],
    atom: Literal["N", "C", "CG"],
    neighbors: List[Atom] = None,
):
    """
    Cap with hydrogen, building based on the upstream or downstream residue

    Parameters
    ----------
    parent: Bio.PDB.Residue
        Residue to cap
    template: Bio.PDB.Residue
        Upstream or downstream residue
    atom: str
        Flag for adding to the 'N' or 'C' or 'CG' (IAS) side of the residue
    neighbors: List[Bio.PDB.Atom], optional
        Neighbor atoms used when an amide hydrogen is missing from ``parent``

    Returns
    -------
    res: Bio.PDB.Residue
        Residue containing added hydrogen
    """
    if neighbors is None:
        neighbors = []
    if template is not None:
        if atom == "N":
            pos = scale_hydrogen(parent["N"], template["C"], 1 / 1.32)
        elif atom == "C":
            pos = scale_hydrogen(parent["C"], template["N"], 1.09 / 1.32)
        elif atom == "CG":
            pos = scale_hydrogen(parent["CG"], template["N"], 1.09 / 1.32)
    else:
        if atom == "N":
            CA = parent["CA"]
            N = parent["N"]
            H = None
            if parent.get_resname() == "PRO":
                # Proline does not have an H atom on N-terminus
                H = parent["CD"]
            elif "H" in parent:
                H = parent["H"]
            else:
                for neighbor in neighbors:
                    if neighbor.element == "H":
                        H = neighbor
                        break
                if H is None:
                    raise KeyError(f"No H atom found for {make_res_key(parent)}")
            bis = get_normalized_vector(N, CA) + get_normalized_vector(N, H)
            bis /= np.linalg.norm(bis)
            pos = N.get_coord() - bis
        elif atom == "C":
            CA = parent["CA"]
            C = parent["C"]
            O = parent["O"]
            bis = get_normalized_vector(C, CA) + get_normalized_vector(C, O)
            bis /= np.linalg.norm(bis)
            pos = C.get_coord() - bis * 1.09
        elif atom == "CG":
            CB = parent["CB"]
            CG = parent["CG"]
            OD1 = parent["OD1"]
            bis = get_normalized_vector(CG, CB) + get_normalized_vector(CG, OD1)
            bis /= np.linalg.norm(bis)
            pos = CG.get_coord() - bis * 1.09

    for name in ["H1", "H2", "H3"]:
        if name not in parent:
            break
    atom = Atom(name, pos, 0, 1, " ", name, None, "H")
    parent.add(atom)
    return atom


def build_heavy(chain, parent, template, atom):
    """
    Cap with ACE/NME, building based on the upstream or downstream residue

    Parameters
    ----------
    chain: Bio.PDB.Chain
        Chain with desired residue
    parent: Bio.PDB.Residue
        Residue to cap
    template: Bio.PDB.Residue
        Upstream or downstream residue
    atom: str
        Flag for adding to the "N" (ACE) or "C" (NME) side of the residue

    Returns
    -------
    res: Bio.PDB.Residue
        Residue containing added group
    """

    pos = {"CH3": template["CA"].get_coord()}
    if template.get_resname() == "GLY":
        pos["HH31"] = template["HA2"].get_coord()
        pos["HH32"] = template["HA3"].get_coord()
    else:
        pos["HH31"] = template["HA"].get_coord()
        pos["HH32"] = scale_hydrogen(template["CA"], template["CB"], 1.09 / 1.54)
    
    if atom == "N":
        pos["C"] = template["C"].get_coord()
        pos["O"] = template["O"].get_coord()
        pos["HH33"] = scale_hydrogen(template["CA"], template["N"], 1.09 / 1.46)
    else:
        pos["N"] = template["N"].get_coord()
        if template.get_resname() == "PRO":
            pos["H"] = scale_hydrogen(template["N"], template["CD"], 1 / 1.46)
        else:
            pos["H"] = template["H"].get_coord()
        pos["HH33"] = scale_hydrogen(template["CA"], template["C"], 1.09 / 1.51)

    adj_id = ("H_" + ("NME" if atom == "N" else "ACE"), template.get_id()[1], " ")
    skip = chain.has_id(adj_id) # skip building methyl if already present in adjacent cap
    if skip:
        chain[adj_id].detach_child("HH33")

    name = "ACE" if atom == "N" else "NME"
    res_id = ("H_" + name, template.get_id()[1], " ")
    res = Residue(res_id, name, " ")
    for k, v in pos.items():
        if skip and k in ["CH3", "HH31", "HH32", "HH33"]:
            continue
        res.add(Atom(k, v, 0, 1, " ", k, None, k[0]))
    chain.add(res)
    return res


def check_atom_valence(
    res: Residue,
    tree: NeighborSearch,
    atom: Literal["N", "C", "CG"],
    cn: int,
    backbone: bool = True,
    same_residue: bool = False,
) -> Tuple[bool, List[Atom]]:
    """Check whether an atom already has sufficient bonded neighbors.

    Uses a 1.8 A distance search to find neighbors. Also checks for
    peptide bond partners (C/CG bonded to N, or N bonded to C/CG) when
    ``backbone`` is True.

    Parameters
    ----------
    res : Bio.PDB.Residue.Residue
        Residue containing the atom.
    tree : Bio.PDB.NeighborSearch
        Spatial search object for the structure.
    atom : str
        Atom name to check (``'N'``, ``'C'``, ``'CG'``, etc.).
    cn : int
        Minimum coordination number indicating the atom is already saturated.
    backbone : bool, optional
        If True, treat peptide-bond partners as satisfying valence.
    same_residue : bool, optional
        If True, only count neighbors that belong to ``res``.

    Returns
    -------
    tuple of (bool, list)
        ``(True, neighbors)`` if the atom already has enough neighbors;
        otherwise ``(False, neighbors)``.
    """
    neighbors = tree.search(res[atom].get_coord(), radius=1.8)
    if same_residue:
        neighbors = [n for n in neighbors if n.get_parent() == res]
    check_flag = False
    if len(neighbors) > cn:
        check_flag = True
    elif backbone:
        for neighbor in neighbors:
            if neighbor.get_name() in ["C", "CG"] and atom == "N":
                check_flag = True
            elif neighbor.get_name() == "N" and atom in ["C", "CG"]:
                check_flag = True
    return check_flag, neighbors


def primed_atom_names(name: str) -> tuple:
    """Return PDB atom-name variants for primed sugar atoms (O3'/O3*)."""
    if name.endswith("'"):
        return (name, name[:-1] + "*")
    if name.endswith("*"):
        return (name, name[:-1] + "'")
    return (name,)


def get_res_atom(res: Residue, name: str) -> Optional[Atom]:
    """Get an atom from ``res``, accepting both ``'`` and ``*`` primed names."""
    for candidate in primed_atom_names(name):
        if res.has_id(candidate):
            return res[candidate]
    return None


def has_res_atom(res: Residue, name: str) -> bool:
    """Return True if ``res`` has ``name`` under either primed spelling."""
    return get_res_atom(res, name) is not None


def o3prime_atom_name(res: Residue) -> Optional[str]:
    """Return the in-residue O3'/O3* atom name, or None if absent."""
    for candidate in primed_atom_names("O3'"):
        if res.has_id(candidate):
            return candidate
    return None


def atom_coordination_is_one(
    res: Residue, tree: NeighborSearch, atom: str
) -> bool:
    """True when ``atom`` has formal coordination number 1 (unsaturated).

    Uses :func:`check_atom_valence` with ``cn=2`` and ``backbone=False``:
    self + one bonded neighbor → unsaturated; an extra H/metal saturates it.
    """
    if not res.has_id(atom):
        return False
    check_flag, _ = check_atom_valence(res, tree, atom, 2, backbone=False)
    return not check_flag


def phosphate_terminal_pair_charge(
    res: Residue, tree: NeighborSearch, o1: str, o2: str
) -> int:
    """Return -1 when both terminal phosphate oxygens have CN == 1."""
    if atom_coordination_is_one(res, tree, o1) and atom_coordination_is_one(
        res, tree, o2
    ):
        return -1
    return 0


def polymer_nucleotide_charge(res: Residue, tree: NeighborSearch) -> int:
    """Formal charge for A/C/G/U or MGT absent from the Protoss ligand SDF."""
    resname = res.get_resname().strip()
    c = 0
    if resname in RNA_POLYMER_RESNAMES:
        c += phosphate_terminal_pair_charge(res, tree, "OP1", "OP2")
    elif resname in MGT_RESNAMES:
        for o1, o2 in (("O1A", "O2A"), ("O1B", "O2B"), ("O1G", "O2G")):
            c += phosphate_terminal_pair_charge(res, tree, o1, o2)
        if res.has_id("N7"):
            c += 1
    return c


def _unit_vec(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-9:
        raise ValueError("zero-length vector")
    return v / n


def _perpendicular(a: np.ndarray) -> np.ndarray:
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, a)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    return _unit_vec(ref - np.dot(ref, a) * a)


def o3prime_has_hydrogen(parent: Residue) -> bool:
    """True if O3' already carries an in-residue hydroxyl hydrogen."""
    o3 = get_res_atom(parent, "O3'")
    if o3 is None:
        return False
    if has_res_atom(parent, "HO3'") or has_res_atom(parent, "H3'"):
        return True
    o = o3.get_coord()
    for atom in parent.get_atoms():
        if atom.element not in ("H", "D"):
            continue
        if float(np.linalg.norm(atom.get_coord() - o)) < 1.2:
            return True
    return False


def build_o3prime_hydrogen(parent: Residue) -> Optional[Atom]:
    """Cap a dangling O3' with hydrogen (3'-OH link-atom cap)."""
    if o3prime_has_hydrogen(parent):
        return None
    o3 = get_res_atom(parent, "O3'")
    c3 = get_res_atom(parent, "C3'")
    if o3 is None or c3 is None:
        return None

    o_coord = o3.get_coord()
    a = _unit_vec(c3.get_coord() - o_coord)
    direction = _unit_vec(_TET_COS * a + _TET_SIN * _perpendicular(a))
    pos = o_coord + O_H_BOND * direction
    name = "HO3'"
    atom = Atom(name, pos, 0, 1, " ", name, None, "H")
    parent.add(atom)
    return atom


def ias_cg_isopeptide_linked(res: Residue, tree: NeighborSearch) -> bool:
    """True if IAS CG is amide-linked to another residue's N (isopeptide)."""
    if not res.has_id("CG"):
        return False
    for neighbor in tree.search(res["CG"].get_coord(), radius=1.8):
        if neighbor.get_parent() is res:
            continue
        if neighbor.element == "N" and neighbor.get_name() == "N":
            return True
    return False


def hetero_residue_formal_charge(
    res: Residue,
    tree: NeighborSearch,
    n_terminals: set,
) -> int:
    """Formal charge for sphere-0 heteros missing from the Protoss SDF.

    Mirrors the AA-loop N-terminus / OXT rules so residues such as IAS keep
    their α-carboxylate (-1) when Protoss omits them from ligands.sdf. IAS
    also gets ASP-like sidechain -1 only when OD1/OD2 are present and CG is
    not isopeptide-linked to another residue.
    """
    c = 0
    res_id = res.get_full_id()
    resname = res.get_resname().strip()
    if res.has_id("N") and res_id in n_terminals and (resname != "PRO" or res.has_id("H")):
        c += 1
    elif res.has_id("N"):
        check_flag, _ = check_atom_valence(res, tree, "N", 4, backbone=False)
        if check_flag:
            c += 1
    if res.has_id("OXT"):
        check_flag, _ = check_atom_valence(res, tree, "OXT", 2, backbone=False)
        if not check_flag:
            c -= 1
    if (
        resname == "IAS"
        and res.has_id("OD1")
        and res.has_id("OD2")
        and all(not res.has_id(h) for h in ["HD2", "HOD1", "HOD2"])
        and not ias_cg_isopeptide_linked(res, tree)
    ):
        c -= 1
    return c


def cap_chains(
    model: Model,
    residues: Set[Residue],
    capping: int,
    RGP_atoms: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    ligand_charge: Optional[dict] = None,
) -> Set[Residue]:
    """
    Cap chain breaks for a set of extracted residues

    Parameters
    ----------
    model: Bio.PDB.Model
        Protein structure model
    residues: set
        Set of residues
    capping: int
        Flag for capping group, H (1) or ACE/NME (2). Non-ligand A/C/G/U
        always receive O3'–H caps when ``capping`` is non-zero.
    RGP_atoms: dict, optional
        RGP atom information used to place hydrogens at missing R# sites
    ligand_charge: dict, optional
        Protoss ligand charge map; residues present here (including oligomer
        members) are skipped for polymer nucleotide O3' capping.

    Returns
    -------
    cap_residues: set
        Set of residues containing added groups
    """
    if RGP_atoms is None:
        RGP_atoms = {}
    if ligand_charge is None:
        ligand_charge = {}
    ligand_keys = ligand_charge.keys()
    orig_chains = {}
    for chain in model:
        orig_chains[chain.get_id()] = chain.get_unpacked_list()

    cap_residues = set()

    cluster_atom_list = []
    for res in residues:
        cluster_atom_list += list(res.get_atoms())
    cluster_tree = NeighborSearch(cluster_atom_list)

    for res in list(sorted(residues)):
        res_key = make_res_key(res)
        if res_key in RGP_atoms:
            for RGP_atom_info in RGP_atoms[res_key].values():
                if RGP_atom_info.get("atom") not in cluster_atom_list:
                    bond_vector = RGP_atom_info["atom"].get_coord() - RGP_atom_info["linking_atom_coord"]
                    norm_bond_vector = bond_vector / np.linalg.norm(bond_vector)
                    linking_element = RGP_atom_info["linking_atom"].element
                    bond_length = HX_BOND_LENGTH.get(linking_element, 1.09)
                    pos = RGP_atom_info["linking_atom_coord"] + norm_bond_vector * bond_length
                    name = "H0"
                    for i in range(100):
                        if f"H{i}" not in res:
                            name = f"H{i}"
                            break
                    res.add(Atom(name, pos, 0, 1, " ", name, None, "H"))

        res_id = res.get_full_id()
        resname = res.get_resname().strip()
        # Polymer RNA bases missing from Protoss ligands: H-cap dangling O3'.
        if resname in RNA_POLYMER_RESNAMES and not residue_in_ligands(
            resname, res_id, False, ligand_keys
        ):
            o3_name = o3prime_atom_name(res)
            if o3_name is not None and atom_coordination_is_one(
                res, cluster_tree, o3_name
            ):
                ho3 = build_o3prime_hydrogen(res)
                if ho3 is not None:
                    cap_residues.add(ho3)

        if not (
            (Polypeptide.is_aa(res) and res.get_id()[0] == " ") # normal amino acid
            or res.get_resname() == "IAS"                       # IAS
        ):
            continue

        chain = model[res_id[2]]
        chain_list = orig_chains[chain.get_id()]
        ind = chain_list.index(res)

        N_capped_flag = False
        if ind > 0:
            pre = chain_list[ind - 1]
            if (
                pre.get_id()[1] == res_id[3][1] - 1
                and pre.get_id()[0] == " "
                and pre not in residues
                and Polypeptide.is_aa(pre)
            ):  # ignores hetero residues
                if capping == 1:
                    cap_residues.add(build_hydrogen(res, pre, "N"))
                else:
                    cap_residues.add(build_heavy(chain, res, pre, "N"))
                N_capped_flag = True
        if not N_capped_flag:
            check_flag, neighbors = check_atom_valence(res, cluster_tree, "N", 3)
            if not check_flag:
                cap_residues.add(build_hydrogen(res, None, "N", neighbors))

        C_capped_flag = False
        if res.get_resname() == "IAS":
            C_name = "CG"
        else:
            C_name = "C"
        if ind < len(chain_list) - 1:
            nxt = chain_list[ind + 1]
            if (
                nxt.get_id()[1] == res_id[3][1] + 1
                and nxt.get_id()[0] == " "
                and nxt not in residues
                and Polypeptide.is_aa(nxt)
            ):
                if capping == 1:
                    cap_residues.add(build_hydrogen(res, nxt, C_name))
                else:
                    cap_residues.add(build_heavy(chain, res, nxt, C_name))
                C_capped_flag = True
        if not C_capped_flag:
            check_flag, neighbors = check_atom_valence(res, cluster_tree, C_name, 3)
            if not check_flag:
                cap_residues.add(build_hydrogen(res, None, C_name, neighbors))

    return cap_residues


def write_pdbs(io, sphere, out):
    """
    Write coordination sphere to PDB file.

    Parameters
    ----------
    io: Bio.PDB.PDBIO
        Bio.PDB writer
    sphere: list
        List of coordination sphere residues
    out: str
        Path to output PDB file
    """

    class ResSelect(Select):
        def accept_residue(self, residue):
            return residue in sphere

    io.save(out, ResSelect())


def residue_in_ligands(resname, resid, res_is_aa, ligand_keys):
    """Check whether a residue matches any key in the ligand charge dictionary.

    For amino acids, an exact key match is required. For non-amino-acid
    residues, the key may be part of a space-separated oligomer key.

    Parameters
    ----------
    resname : str
        Three-letter residue name.
    resid : tuple
        Residue full ID tuple from BioPython.
    res_is_aa : bool
        Whether the residue is an amino acid.
    ligand_keys : iterable of str
        Ligand charge dictionary keys.

    Returns
    -------
    bool
        True if the residue matches a ligand key.
    """
    res_key = f"{resname}_{resid[2]}{resid[3][1]}"
    if res_is_aa:
        return res_key in ligand_keys
    else:
        for ligand_key in ligand_keys:
            ligand_res_keys = ligand_key.split()
            if res_key in ligand_res_keys:
                return True
        return False


def check_disulfide(res: Residue, tree: NeighborSearch):
    """Detect whether a cysteine residue is involved in a disulfide bond.

    A disulfide is detected when both the query CYS and a neighboring CYS
    lack an HG atom and their SG atoms are within 2.5 A.

    Parameters
    ----------
    res : Bio.PDB.Residue.Residue
        A cysteine residue to check.
    tree : Bio.PDB.NeighborSearch
        Spatial search object for the structure.

    Returns
    -------
    bool
        True if the residue is part of a disulfide bond.
    """
    if not res.has_id("HG"):
        SG = res["SG"]
        neighbors = tree.search(SG.get_coord(), radius=2.5)
        for neighbor in neighbors:
            if neighbor != SG and neighbor.get_name() == "SG":
                neighbor_res = neighbor.get_parent()
                if neighbor_res.get_resname() == "CYS" and not neighbor_res.has_id("HG"):
                    return True
    return False


def compute_charge(
    spheres,
    structure,
    ligand_charge,
    center_residue,
    residues=None,
    RGP_atoms=None,
):
    """
    Computes the total charge of coordinating AAs

    Parameters
    ----------
    spheres: list of sets
        Sets of residues separated by spheres
    structure: Bio.PDB.Structure
        The protein structure
    ligand_charge: dict
        Key, value pairs of ligand names and charges
    center_residue: CenterResidue
        The residues to use as the cluster center
    residues: set, optional
        All residues in the cluster; defaults to the union of ``spheres``
    RGP_atoms: dict, optional
        RGP atom information used to avoid double-counting CYS charges

    Returns
    -------
    charge: list
        Total charge of AAs in each sphere
    """
    if RGP_atoms is None:
        RGP_atoms = {}
    if residues is None:
        residues = reduce(lambda x, y: x | y, spheres) if spheres else set()

    # Identifying N-terminal and C-terminal residues for each chain
    n_terminals = set()
    c_terminals = set()
    # Loop over the residues to get first and last as indices may be different
    for chain in structure.get_chains():
        chain_residues = list(chain.get_residues())
        if chain_residues:
            n_terminals.add(chain_residues[0].get_full_id())
            c_terminals.add(chain_residues[-1].get_full_id())

    pos = {
        "ARG": ["HE", "HH11", "HH12", "HH21", "HH22"],
        "LYS": ["HZ1", "HZ2", "HZ3"],
        "HIS": ["HD1", "HD2", "HE1", "HE2"],
        "HIP": ["HD1", "HD2", "HE1", "HE2"],
        "HID": ["HD1", "HD2", "HE1", "HE2"],
        "MLZ": [],
        "M3L": []
    }
    neg = {
        "ASP": ["HD2", "HOD1", "HOD2"],
        "GLU": ["HE2", "HOE1", "HOE2"],
        "CYS": ["HG"],
        "TYR": ["HH"],
        "OCS": [],
        "CSD": ["HD1", "HD2"],
        "KCX": ["HQ1", "HQ2", "HOQ1", "HOQ2"],
        "HIS": ["HD1", "HE2"]
    }

    charge = []
    start_sphere_id = 0 if center_residue.mode == "strict" else 1
    if start_sphere_id == 1:
        charge.append(0)
    
    cluster_atom_list = []
    for s in spheres:
        for res in s:
            cluster_atom_list.extend(list(res.get_atoms()))
    cluster_tree = NeighborSearch(cluster_atom_list)
    res_keys = set(make_res_key(res) for res in residues)

    # Sphere-0 hetero residues missing from the Protoss SDF still need a charge
    # entry so oligomer / center accounting stays consistent downstream. Amino
    # acids remain in the AA charge loop below so CenterResidue strict/extended
    # semantics from hjkgrp main are preserved.
    if spheres:
        s0 = spheres[0]
        sphere_tree = NeighborSearch([atom for res in s0 for atom in res.get_atoms()])
        for res in s0:
            res_id = res.get_full_id()
            resname = res.get_resname()
            res_is_aa = Polypeptide.is_aa(res)
            if res_is_aa:
                continue
            if not residue_in_ligands(resname, res_id, res_is_aa, ligand_charge.keys()):
                resname_key = resname.strip()
                # Polymer nucleotides belong in the sphere charge loop below,
                # not the ligand CSV map (Protoss already covers true ligands).
                if resname_key in RNA_POLYMER_RESNAMES or resname_key in MGT_RESNAMES:
                    continue
                ligand_charge[make_res_key(res)] = hetero_residue_formal_charge(
                    res, sphere_tree, n_terminals
                )

    for s in spheres[start_sphere_id:]:
        sphere_tree = NeighborSearch([atom for res in s for atom in res.get_atoms()])
        c = 0
        for res in s:
            res_id = res.get_full_id()
            resname = res.get_resname()
            res_is_aa = Polypeptide.is_aa(res)
            if not residue_in_ligands(resname, res_id, res_is_aa, ligand_charge.keys()):
                resname_key = resname.strip()
                if resname_key in RNA_POLYMER_RESNAMES or resname_key in MGT_RESNAMES:
                    delta = polymer_nucleotide_charge(res, sphere_tree)
                    if delta:
                        charge_debug(f"polymer nucleotide {delta}", res)
                    c += delta
                    continue
                # Keep main's `pos and all(H)` gate so residues also listed in
                # `neg` (e.g. deprotonated HIS) can fall through to the neg branch.
                if resname in pos and all(res.has_id(h) for h in pos[resname]):
                    charge_debug("pos res +1", res)
                    c += 1
                elif resname == "LYS":
                    check_flag, _ = check_atom_valence(
                        res, sphere_tree, "NZ", 4, backbone=False, same_residue=True
                    )
                    if check_flag:
                        charge_debug("LYS +1", res)
                        c += 1
                elif resname in neg and all(not res.has_id(h) for h in neg[resname]):
                    RGP_flag = False
                    if resname == "CYS":
                        if check_disulfide(res, cluster_tree):
                            RGP_flag = True
                        elif "SG" in res:
                            for name, RGP_atom_list in RGP_atoms.items():
                                for RGP_atom_info in RGP_atom_list.values():
                                    RGP_atom = RGP_atom_info.get("atom")
                                    if res["SG"] == RGP_atom and name in res_keys:
                                        RGP_flag = True
                                        break
                                if RGP_flag:
                                    break
                    if not RGP_flag:
                        charge_debug("neg res -1", res)
                        c -= 1
                if res_is_aa and resname != "PRO" and all(not res.has_id(h) for h in ["H", "H2"]):
                    # Deprotonated backbone amide (missing H/H2) is formally -1.
                    # Do not gate on check_atom_valence: peptide C/CA neighbors are
                    # expected, and metal–N coordination must not cancel this charge.
                    charge_debug("backbone N -1", res)
                    c -= 1

                # Check for charged N-terminus (NH3+ / Pro NH2+).
                # Protoss may leave a "fake" N-terminus with a single amide-like H;
                # H-capping then yields neutral NH2. NeighborSearch counts the N
                # itself, so CN > 4 marks NH3+ / Pro-NH2+ (5 neighbors) but not
                # capped NH2 (4 neighbors).
                if res_id in n_terminals and res.has_id("N"):  # exclude sugar chain terminus
                    check_flag, _ = check_atom_valence(
                        res, cluster_tree, "N", 4, backbone=False
                    )
                    if check_flag:
                        charge_debug("sphere 1+ N terminal +1", res)
                        c += 1

                # Check for charged C-terminus
                if res.has_id("OXT"):
                    check_flag, _ = check_atom_valence(
                        res, sphere_tree, "OXT", 2, backbone=False
                    )
                    if not check_flag:
                        charge_debug("sphere 1+ C terminal -1", res)
                        c -= 1

        charge.append(c)
    return charge


def count_residues(spheres):
    """
    Counts the frequency of coordinating residues

    Parameters
    ----------
    spheres: list of sets
        Sets of residues separated by spheres

    Returns
    -------
    count: list of dicts
        Frequency table by sphere
    """
    count = []
    for s in spheres[1:]:
        c = {}
        for res in s:
            c[res.get_resname()] = c.get(res.get_resname(), 0) + 1
        count.append(c)
    return count


def make_res_key(res):
    """Format a residue as a ``'RESNAME_CHAINID'`` string key (e.g., ``'FE_A199'``)."""
    resname = res.get_resname()
    resid = res.get_id()[1]
    chainid = res.get_parent().get_id()
    return f"{resname}_{chainid}{resid}"    


def complete_oligomer(ligand_keys, model, residues, spheres, include_ligands) -> List[Residue]:
    """Ensure that partially included oligomeric ligands are fully added.

    If any residue of a multi-residue ligand (oligomer) is present in the
    extracted spheres, all remaining residues of that oligomer are added to
    avoid unpredictable charge errors.

    Parameters
    ----------
    ligand_keys : iterable of str
        Ligand charge dictionary keys (space-separated for oligomers).
    model : Bio.PDB.Model.Model
        Full protein structure model.
    residues : set
        Current set of extracted residues (modified in place).
    spheres : list of set
        Sphere-separated residue sets (modified in place).
    include_ligands : int
        Ligand inclusion level (0 = first sphere only, 1 = non-water, 2 = all,
        3 = center/standard AA/water).

    Returns
    -------
    list of Bio.PDB.Residue.Residue
        Oligomer monomers that should be protected from pruning.
    """
    ligand_res_found = dict()
    oligomer_found = dict()
    for ligand_key in ligand_keys:
        ligand_res_keys = ligand_key.split()
        if len(ligand_res_keys) == 1:
            continue
        oligomer_found[ligand_key] = dict()
        for ligand_res_key in ligand_res_keys:
            ligand_res_found[ligand_res_key] = {
                "sphere": -1,
                "oligomer": ligand_key
            }
            oligomer_found[ligand_key][ligand_res_key] = False
    if not oligomer_found:
        return []
    kept_monomers = []
    for i, sphere in enumerate(spheres):
        if include_ligands == 0 and i > 0:
            break
        for res in sphere:
            res_key = make_res_key(res)
            if res_key in ligand_res_found and not Polypeptide.is_aa(res):
                ligand_res_found[res_key]["sphere"] = i
                oligomer = ligand_res_found[res_key]["oligomer"]
                oligomer_found[oligomer][res_key] = True
    for chain in model:
        for res in chain.get_unpacked_list():
            res_key = make_res_key(res)
            if (
                res_key in ligand_res_found and 
                not Polypeptide.is_aa(res)
            ):
                oligomer = ligand_res_found[res_key]["oligomer"]
                found_sphere = ligand_res_found[res_key]["sphere"]
                if any(oligomer_found[oligomer].values()):
                    kept_monomers.append(res)
                    if found_sphere < 0:
                        if include_ligands == 0:
                            spheres[0].add(res)
                        else:
                            spheres[-1].add(res)
                        residues.add(res)
                        print(f"To avoid unpredictable charge error, {res_key} in {oligomer} is added to spheres")
    return kept_monomers


def find_RGP_atoms(structure: Structure, RGP_atoms: Dict[str, Dict[int, Dict[str, Any]]]) -> None:
    """Match SDF RGP / linking-atom coordinates onto PDB atoms in place."""
    if not RGP_atoms:
        return
    for atom in structure.get_atoms():
        atom_coord = atom.get_coord()
        for RGP_atom_list in RGP_atoms.values():
            for RGP_atom_info in RGP_atom_list.values():
                if "coord" in RGP_atom_info and np.allclose(
                    atom_coord, RGP_atom_info["coord"], atol=1e-3
                ):
                    RGP_atom_info["atom"] = atom
                if "linking_atom_coord" in RGP_atom_info and np.allclose(
                    atom_coord, RGP_atom_info["linking_atom_coord"], atol=1e-3
                ):
                    RGP_atom_info["linking_atom"] = atom


def extract_clusters(
    path,
    out,
    center_residue: CenterResidue,
    sphere_count=2,
    first_sphere_radius=4.0,
    max_atom_count=None,
    merge_cutoff=0.0,
    smooth_method="box_plot",
    ligands=[],
    capping=1,
    charge=True,
    ligand_charge=dict(),
    count=True,
    xyz=True,
    hetero_pdb=False,
    include_ligands=2,
    cluster_name_template=None,
    RGP_atoms=None,
    **smooth_params
):
    """Extract active site coordination spheres using Voronoi tessellation.

    The main entry point for cluster extraction. Identifies center residues,
    builds coordination spheres using Voronoi neighbors, applies optional
    capping, and writes output files (PDB, XYZ, charge.csv, count.csv).

    Parameters
    ----------
    path : str
        Path to the input PDB file.
    out : str
        Path to the output directory.
    center_residue : CenterResidue
        Definition of residues to use as the cluster center.
    sphere_count : int, optional
        Number of coordination spheres to extract (default 2).
    first_sphere_radius : float, optional
        Distance cutoff in angstroms for the first sphere (default 4.0).
    max_atom_count : int, optional
        Maximum atom count; residues are pruned if exceeded (default None).
    merge_cutoff : float, optional
        Distance cutoff for merging nearby centers (default 0.0).
    smooth_method : str, optional
        Sphere smoothing method: ``'box_plot'``, ``'dbscan'``, or
        ``'dummy_atom'`` (default ``'box_plot'``).
    ligands : list, optional
        Additional ligand residue names to include (default []).
    capping : int, optional
        Capping mode: 0 = none, 1 = hydrogen, 2 = ACE/NME (default 1).
    charge : bool, optional
        If True, write amino acid charges to ``charge.csv`` (default True).
    ligand_charge : dict, optional
        Mapping of ligand IDs to formal charges (default {}).
    count : bool, optional
        If True, write residue counts to ``count.csv`` (default True).
    xyz : bool, optional
        If True, write XYZ coordinate files (default True).
    hetero_pdb : bool, optional
        If True, include HETATM records in combined PDB (default False).
    include_ligands : int, optional
        Ligand inclusion mode: 0 = first sphere only unless in ``ligands``,
        1 = all non-water, 2 = all (default 2), 3 = center / standard amino
        acids / waters only.
    cluster_name_template : str, optional
        Python format-string controlling cluster directory/file names.
        Defaults to ``None``, which preserves the original behavior of
        naming each cluster after its chain + residue number(s), e.g.
        ``A199``, or ``A1_A2_A3_A4`` for a cluster merging several centers
        (see ``merge_cutoff``). This can produce filenames long enough to
        exceed OS path length limits (notably on Windows) when many
        residues are merged into one center, so a template can be supplied
        instead. The template is evaluated once per cluster with these
        fields available:

        - ``radius``: ``first_sphere_radius``, formatted without a
          trailing ``.0`` (e.g. ``4`` or ``3.5``).
        - ``metal_id``: the original chain + residue string, e.g. ``A199``.
        - ``index``: 1-based count of the cluster within this call
          (1, 2, 3, ...).
        - ``pdb``: the base name of ``out`` (typically the PDB ID).

        For example ``"A_{radius}"`` names clusters like ``A_4``, and
        ``"cluster_{index}"`` names them ``cluster_1``, ``cluster_2``, etc.
        If two clusters in the same call resolve to the same name (e.g.
        several distinct centers sharing one radius), a numeric suffix is
        appended automatically so nothing gets overwritten (``A_4``,
        ``A_4_1``, ``A_4_2``, ...). ``charge.csv`` and ``count.csv`` rows
        are keyed by the final cluster name (not ``metal_id``), since
        downstream QM job creation (``qp create``/``qp submit``) looks up
        each cluster's charge by matching its directory name against
        these files. When a template renames a cluster away from its
        residue-based ``metal_id``, that original identity is instead
        recorded in ``cluster_name_map.csv`` (``cluster_name,metal_id``)
        so it isn't lost.
    RGP_atoms : dict, optional
        Mapping of ligand keys to RGP atom metadata from
        :func:`qp.protonate.ligand_prop.collect_RGP_atoms`. Used for RGP
        hydrogen capping and CYS charge corrections.
    **smooth_params
        Additional parameters for the smoothing method.

    Returns
    -------
    list of str
        Paths to the generated cluster directories (e.g.,
        ``['out/A199', 'out/B350']``, or ``['out/A_4', 'out/A_4_1']`` with
        ``cluster_name_template="A_{radius}"``).
    """
    if RGP_atoms is None:
        RGP_atoms = {}
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", path)
    io = PDBIO()
    io.set_structure(structure)

    model = structure[0]
    neighbors = voronoi(model, center_residue, ligands, smooth_method, out,**smooth_params)

    centers = get_center_residues(model, center_residue, merge_cutoff)

    aa_charge = {}
    res_count = {}
    cluster_paths = []
    cluster_names_used = {}
    cluster_name_map = {}
    pdb_name = os.path.basename(os.path.normpath(out))
    for index, c in enumerate(centers, start=1):
        metal_id, residues, spheres = get_next_neighbors(
            c, neighbors, sphere_count, ligands, first_sphere_radius, smooth_method, include_ligands, **smooth_params
        )
        kept_monomers = complete_oligomer(ligand_charge, model, residues, spheres, include_ligands)

        if cluster_name_template:
            name_fields = {
                "radius": f"{first_sphere_radius:g}",
                "metal_id": metal_id,
                "index": index,
                "pdb": pdb_name,
            }
            try:
                cluster_name = cluster_name_template.format(**name_fields)
            except (KeyError, IndexError) as e:
                raise ValueError(
                    f"Invalid cluster_name_template {cluster_name_template!r}: "
                    f"unknown field {e}. Available fields: {sorted(name_fields)}"
                ) from e
        else:
            cluster_name = metal_id

        # Guard against two clusters resolving to the same name (e.g. several
        # distinct centers sharing the same radius) so nothing overwrites.
        if cluster_name in cluster_names_used:
            cluster_names_used[cluster_name] += 1
            cluster_name = f"{cluster_name}_{cluster_names_used[cluster_name]}"
        else:
            cluster_names_used[cluster_name] = 0

        cluster_path = f"{out}/{cluster_name}"
        cluster_paths.append(cluster_path)
        os.makedirs(cluster_path, exist_ok=True)

        if cluster_name != metal_id:
            cluster_name_map[cluster_name] = metal_id

        find_RGP_atoms(structure, RGP_atoms)
        if max_atom_count is not None:
            prune_atoms(c, residues, spheres, max_atom_count, ligands, kept_monomers)
        if count:
            res_count[cluster_name] = count_residues(spheres)
        if capping:
            cap_residues = cap_chains(
                model, residues, capping, RGP_atoms, ligand_charge=ligand_charge
            )
            if capping == 2:
                spheres[-1] |= cap_residues
        if charge:
            aa_charge[cluster_name] = compute_charge(
                spheres,
                structure,
                ligand_charge,
                center_residue,
                residues=residues,
                RGP_atoms=RGP_atoms,
            )

        sphere_paths = []
        for i, s in enumerate(spheres):
            sphere_path = f"{cluster_path}/{i}.pdb"
            sphere_paths.append(sphere_path)
            write_pdbs(io, s, sphere_path)
        if capping:
            for cap in cap_residues:
                cap.get_parent().detach_child(cap.get_id())
        if xyz:
            struct_to_file.to_xyz(f"{cluster_path}/{cluster_name}.xyz", *sphere_paths)
            struct_to_file.combine_pdbs(f"{cluster_path}/{cluster_name}.pdb", center_residue, *sphere_paths, hetero_pdb=hetero_pdb)

    if charge:
        with open(f"{out}/charge.csv", "w") as f:
            f.write(f"Name,{','.join(str(i) for i in range(sphere_count + 1))}\n")
            for k, v in sorted(aa_charge.items()):
                f.write(k)
                for s in v:
                    f.write(f",{s}")
                f.write(f"\n")

    if count:
        with open(f"{out}/count.csv", "w") as f:
            f.write(f"Name,{','.join(str(i + 1) for i in range(sphere_count))}\n")
            for k, v in sorted(res_count.items()):
                f.write(k)
                for sphere in v:
                    s = ", ".join(f"{r} {c}" for r, c in sorted(sphere.items()))
                    f.write(f',"{s}"')
                f.write("\n")

    if cluster_name_map:
        # Only written when cluster_name_template renamed at least one
        # cluster away from its residue-based metal_id, so the mapping
        # back to the original chain/residue identity isn't lost.
        with open(f"{out}/cluster_name_map.csv", "w") as f:
            f.write("cluster_name,metal_id\n")
            for cluster_name, metal_id in sorted(cluster_name_map.items()):
                f.write(f"{cluster_name},{metal_id}\n")

    return cluster_paths






import numpy as np
from scipy.spatial import Voronoi
import matplotlib.pyplot as plt

def format_plot() -> None:
    """General plotting parameters for the Kulik Lab."""
    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams["xtick.major.pad"] = 5
    plt.rcParams["ytick.major.pad"] = 5
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["svg.fonttype"] = "none"

def plot_voronoi_2d(points, points_count, output_path, x_threshold=(22, 32), y_threshold=(15, 25)):
    """
    Plot a 2D Voronoi diagram from the atomic coordinates, excluding unbounded vertices and filtering by x and y thresholds.
    
    Parameters
    ----------
    points: list or numpy.array
        The 3D coordinates of the points (atoms) to be tessellated.
        Only the xy-plane projection will be plotted.
    points_count: int
        The number of real atoms (excluding dummy atoms).
    output_path: str
        The directory where the plot will be saved.
    x_threshold: tuple of floats
        The min and max thresholds for the x-coordinate filtering.
    y_threshold: tuple of floats
        The min and max thresholds for the y-coordinate filtering.

    Notes
    -----
    Zoom-Out ->  x_threshold=(-22, 83), y_threshold=(-20, 62)
    Zoom-In  ->  x_threshold=(22, 32), y_threshold=(15, 25)

    """
    # Ensure points is a NumPy array
    points = np.array(points)
    
    # Compute the Voronoi tessellation using only xy coordinates
    vor = Voronoi(points[:, :2])
    
    # Plot the Voronoi diagram
    format_plot()
    fig, ax = plt.subplots()

    # Filter vertices to exclude extreme values for both x and y coordinates
    x_min, x_max = x_threshold
    y_min, y_max = y_threshold

    valid_vertices = vor.vertices[
        (vor.vertices[:, 0] >= x_min) & (vor.vertices[:, 0] <= x_max) & 
        (vor.vertices[:, 1] >= y_min) & (vor.vertices[:, 1] <= y_max)
    ]
    
    ax.plot(valid_vertices[:, 0], valid_vertices[:, 1], 'o', markersize=3, color='red', zorder=1)

    # Plot ridges (edges) of the Voronoi cells, but only for real atoms
    for ridge in vor.ridge_vertices:
        if all(v >= 0 for v in ridge):  # Ignore unbounded vertices
            ridge_vertices = vor.vertices[ridge]
            # Filter the ridge by the x and y thresholds
            if (np.all(ridge_vertices[:, 0] >= x_min) and np.all(ridge_vertices[:, 0] <= x_max) and 
                np.all(ridge_vertices[:, 1] >= y_min) and np.all(ridge_vertices[:, 1] <= y_max)):
                ax.plot(ridge_vertices[:, 0], ridge_vertices[:, 1], 'k-')

    # Plot the original points (real atoms)
    # ax.plot(points[:points_count, 0], points[:points_count, 1], 'o', markersize=3, color='red', label='Atoms')

    ax.set_xlabel('X Coordinate (Å)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Y Coordinate (Å)', fontsize=10, fontweight='bold')
    ax.set_aspect('equal')

    # Save the plot to the output path
    output_file = os.path.join(output_path, '2D_filtered_voronoi.png')
    plt.savefig(output_file, bbox_inches="tight", dpi=600)
    plt.close()

