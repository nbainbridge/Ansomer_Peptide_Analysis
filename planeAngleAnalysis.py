"""
PlaneAngleWesternEastern.py
===============================
Author  : Nicholas Bainbridge (nbainbr@student.ubc.ca)
Requires: Python 3, NumPy
 
Description
-----------
Measures how a bond vector (C6->N7) is tilted relative to two different
ring planes, on the "Western" and "Eastern" sides of the ring, across a set
of MD snapshot PDB files.
 
    Western plane : atoms CA2, CA7, C6     Vector : C6 -> N7
    Eastern plane : atoms CA2, CA1, C6     Vector : C6 -> N7
 
For each plane, the plane normal is found from the three atoms, and the
oriented angle between the C6->N7 vector and the plane is returned in
degrees (0 to 180):
 
    0 - 90°   : vector points to the same side as the plane normal
                (90° = along the normal, 0° = lying in the plane)
    90 - 180° : vector points to the opposite side of the plane
 
Which side counts as "positive" depends on the order of the three plane
atoms, so keep that order fixed when comparing results.
 
The script walks every subfolder of `base_dir` (e.g. one subfolder per
simulation/condition) and writes a Western and an Eastern results file
into each subfolder.
 
Sample input
------------
    torsionAngleAnalysis/        <- base_dir
        Water/                   <- one subfolder per condition
            frame0.pdb
            frame500.pdb
            ...
        DMSO/
            frame0.pdb
            ...
 
Each PDB must contain atoms named CA1, CA2, CA7, C6 and N7.
 
Sample output
-------------
Two files written into each subfolder:
 
    PlaneAngleAnalysis/Water/WesternSide_data_2026-09-29_14-05-12.txt
    PlaneAngleAnalysis/Water/EasternSide_data_2026-09-29_14-05-12.txt
 
Contents (PDB file name, oriented vector-plane angle in degrees):
 
    frame0 34.918271
    frame500 41.203385
    frame1000 128.550914
    ...
 
Console:
 
    Processed Water
    Processed DMSO
"""
 
import numpy as np
from pathlib import Path
from datetime import datetime
import os
 
# ---- USER SETTING: parent folder containing one subfolder per condition ----
base_dir = Path(r"PATH")
    
def get_plane(planePointOne, planePointTwo, planePointThree):
    """
    Return the normal vector of the plane through three points.
 
    The normal is (P2 - P1) x (P3 - P1), so its direction depends on the
    order of the points. It is NOT normalized here.
 
    Raises
    ------
    ValueError
        If the three points are collinear (no unique plane).
    """
    planePointOne = np.asarray(planePointOne, dtype=float)
    planePointTwo = np.asarray(planePointTwo, dtype=float)
    planePointThree = np.asarray(planePointThree, dtype=float)
 
    # Two in-plane vectors that share P1 as their origin
    plane_vector_one = planePointTwo - planePointOne
    plane_vector_two= planePointThree - planePointOne
 
    # Their cross product is perpendicular to the plane
    normal = np.cross(plane_vector_one, plane_vector_two)
    if np.linalg.norm(normal) == 0:
        raise ValueError("Plane points are collinear.")
    return normal
 
 
def get_vector(vectorPointOne, vectorPointTwo):
    """
    Return the vector from point one to point two (P2 - P1).
 
    Raises
    ------
    ValueError
        If the two points are identical (zero-length vector).
    """
    vectorPointOne= np.asarray(vectorPointOne, dtype=float)
    vectorPointTwo= np.asarray(vectorPointTwo, dtype=float)
 
    vector = vectorPointTwo - vectorPointOne
    if np.linalg.norm(vector) == 0:
        raise ValueError("Vector has zero length.")
    return vector
 
 
def vector_plane_angle(planeNormal, vector, degrees=True):
    """
    Returns an oriented angle between the vector and the plane.
 
    0-90°  : vector points toward the positive plane normal
    90-180°: vector points toward the opposite side of the plane
 
    The orientation depends on the ordering of the three plane atoms.
 
    Parameters
    ----------
    planeNormal : array-like, shape (3,)
        Plane normal, e.g. from get_plane().
    vector : array-like, shape (3,)
        Vector to compare, e.g. from get_vector().
    degrees : bool
        Return degrees (True, default) or radians (False).
    """
    # Work with unit vectors
    planeNormal = planeNormal / np.linalg.norm(planeNormal)
    vector = vector / np.linalg.norm(vector)
 
    # signed cosine with the plane normal
    cos_theta = np.clip(np.dot(vector, planeNormal), -1.0, 1.0)
 
    # angle between vector and normal (0-180°)
    theta_normal = np.degrees(np.arccos(cos_theta))
 
    # convert to oriented angle with plane
    if theta_normal <= 90:
        angle = 90 - theta_normal
    else:
        angle = 90 + (theta_normal - 90)
 
    return angle if degrees else np.radians(angle)
 
 
def load_pdb(filename):
    """
    Read a PDB file into a dictionary of {atom_name: (x, y, z)}.
 
    Only ATOM/HETATM lines are read. If an atom name appears more than once,
    the last occurrence overwrites the earlier ones.
    """
    atoms = {}
 
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # PDB is a fixed-width format; these column slices are standard
                name=line[12:16].strip()
                atoms[name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
 
    return atoms
 
 
def write_results(results, output_folder, prefix):
    """
    Write (name, value) pairs to a timestamped text file.
 
    Output file: <output_folder>/<prefix>_data_<YYYY-mm-dd_HH-MM-SS>.txt
    One line per result: "<name> <value to 6 d.p.>"
    """
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
 
    output_file = Path(output_folder) / f"{prefix}_data_{timestamp}.txt"
 
    with open(output_file, "w") as f:
        for name, angle in results:
            f.write(f"{name} {angle:.6f} \n")
 
 
if __name__ == "__main__":
 
    os.system("cls")  # clear the console (Windows only)

    # Loop over each subfolder in base_dir
    for subdir in base_dir.iterdir():
        if not subdir.is_dir():
            continue  # skip files in base directory
    
        results_WesternSide = []
        results_EasternSide = []
    
        # Find PDB files inside this subfolder
        pdb_files = list(subdir.glob("*.pdb"))
        if not pdb_files:
            continue  # skip folders with no pdbs
    
        for pdb_file in pdb_files:
            atoms = load_pdb(pdb_file)
    
            #"Western Side" : plane CA2, CA7, C6. Vector defined by: C6 -> N7
            angle_WesternSide = vector_plane_angle(get_plane(atoms["CA2"],atoms["CA7"],atoms["C6"]), get_vector(atoms["C6"],atoms["N7"]))
            results_WesternSide.append((pdb_file.stem, angle_WesternSide))
    
            #"Eastern Side" : plane CA2, CA1, C6. Vector defined by: C6 -> N7
            angle_EasternSide = vector_plane_angle(get_plane(atoms["CA2"],atoms["CA1"],atoms["C6"]), get_vector(atoms["C6"],atoms["N7"]))
            results_EasternSide.append((pdb_file.stem, angle_EasternSide))
    
            # NOTE: these are inside the per-file loop, so the output files are
            # re-written after every PDB (see review notes).
            write_results(results_WesternSide, subdir, "WesternSide")
            write_results(results_EasternSide, subdir, "EasternSide")
            print(f"Processed {subdir.name}")
    
    
        # Quick sanity check of vector_plane_angle with simple coordinates:
        # p1 = (0, 0, 0)
        # p2 = (1, 0, 0)
        # p3 = (0, 1, 0)
        # p4 = (2, 5, 3)
        # print(vector_plane_angle(get_plane(p1,p2,p3), get_vector(p3,p4)))