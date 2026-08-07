"""AMBER ff14SB force field residue name conversion to canonical PDB"""
# residue list from amino12.lib ($AMBERHOME/dat/leap/lib/)

def get_amber_to_pdbcanonical_name():
    """Return a flat {amber_name : canonical_pdb_name} dictionary."""
    return {v:k for k, variants in get_resname_dict().items() for v in variants}

def get_resname_dict():
    """Dictionary of AMBER ff14SB residue names that correspond to each canonical PBD residue name."""
    resname_dict = {
        "HIS" : ["HIE", "HID", "HIP"],
        "ASP" : ["ASH"],
        "CYS" : ["CYM", "CYX"],
        "GLU" : ["GLH"],
        "LYS" : ["LYN"],
    }
    return resname_dict