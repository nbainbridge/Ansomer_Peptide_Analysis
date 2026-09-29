"""
RingAngle.py
============
Author  : Nicholas Bainbridge (nbainbr@student.ubc.ca)
Requires: Python 3, NumPy
 
Description
-----------
Batch-calculates a three-atom bond angle inside the ring across a set of
MD snapshot PDB files.
 
The angle measured is CA7-CA6-CA5, with CA6 as the vertex. It is the
angle between the vectors CA6->CA7 and CA6->CA5, returned in degrees
(0 to 180).
 
The script walks every subfolder of `base_dir` (e.g. one subfolder per
simulation/condition), measures the angle for every PDB in that
subfolder, and writes one results file into that subfolder.
 
Atoms are looked up by their PDB atom name (columns 13-16), so atom names
must be unique within each file.
 
Sample input
------------
    newRingAngle_II/             <- base_dir
        Water/                   <- one subfolder per condition
            frame0.pdb
            frame500.pdb
            ...
        DMSO/
            frame0.pdb
            ...
 
Each PDB must contain atoms named CA5, CA6 and CA7.
 
Sample output
-------------
One file written into each subfolder:
 
    newRingAngle_II/Water/ringAngle_data_2026-09-29_14-05-12.txt
 
Contents (PDB file name, CA7-CA6-CA5 angle in degrees):
    frame0 118.402731
    frame500 121.884209
    frame1000 116.957340
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
# ---- USER SETTING: Atoms for representing the angles NOTE: also need to edit the key in ln 132-134 ----
required_atoms = ["CA6", "CA7", "CA5"]
 
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
 
 
def load_pdb(filename):
    """
    Read a PDB file into a dictionary of {atom_name: np.array([x, y, z])}.
 
    Only ATOM/HETATM lines are read. If an atom name appears more than once,
    the last occurrence overwrites the earlier ones.
    """
    atoms = {}
 
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # PDB is a fixed-width format; these column slices are standard
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
 
                atoms[name] = np.array([x, y, z], dtype=float)
 
    return atoms
 
 
def ringAngle(atoms):
    """
    Calculate the angle CA7-CA6-CA5 (degrees), using CA6 as the vertex.
 
    Parameters
    ----------
    atoms : dict
        {atom_name: coordinates} as returned by load_pdb().
 
    Returns
    -------
    float
        Angle in degrees, 0 to 180.
 
    Raises
    ------
    KeyError
        If CA5, CA6 or CA7 is missing from the PDB.
    """
 

 
    for atom in required_atoms:
        if atom not in atoms:
            raise KeyError(f"Missing atom {atom} in PDB")
 
    p1  = atoms["CA7"]
    vertex = atoms["CA6"]  # vertex
    p2 = atoms["CA5"]
 
    # Vectors from the vertex out to the two other atoms
    v1 = p1 - vertex
    v2 = p2 - vertex
 
    # cos(theta) = (v1 . v2) / (|v1| |v2|)
    dot = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
 
    cos_theta = dot / (norm_v1 * norm_v2)
    # Clip to [-1, 1] so floating-point round-off can't make arccos return NaN
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
 
    angle_rad = np.arccos(cos_theta)
    angle_deg = np.degrees(angle_rad)
 
    return angle_deg
 
 
if __name__ == "__main__":
    os.system("cls")  # clear the console (Windows only)
 
    # Loop over each subfolder in base_dir
    for subdir in base_dir.iterdir():
        if not subdir.is_dir():
            continue  # skip files in base directory
 
        results = []
 
        # Find PDB files inside this subfolder
        pdb_files = list(subdir.glob("*.pdb"))
        if not pdb_files:
            continue  # skip folders with no pdbs
 
        for pdb_file in pdb_files:
            atoms = load_pdb(pdb_file)
            angle = ringAngle(atoms)
            results.append((pdb_file.stem, angle))
 
        # Write results inside the SAME subfolder
        write_results(results, subdir, "ringAngle")
 
        print(f"Processed {subdir.name}")