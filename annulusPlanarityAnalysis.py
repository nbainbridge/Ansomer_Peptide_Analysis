"""
Annulus.py
==========
Author  : Nicholas Bainbridge (nbainbr@student.ubc.ca)
Requires: Python 3, NumPy

Description
-----------
Measures how planar (flat) a ring is across a set of MD snapshots.

For every .pdb file in a folder, all ATOM/HETATM coordinates are read, a
best-fit plane is fitted through them (via SVD), and the root-mean-square
(RMS) deviation of the atoms from that plane is reported. A value of 0 Å
means the atoms are perfectly coplanar; larger values mean more puckering
or bending out of the plane.

Each PDB is expected to contain ONLY the ring atoms of interest (e.g. PDBs
that were split/trimmed down to the macrocycle atoms beforehand). Every
atom in the file is used in the fit.

Sample input
------------
A folder of PDB files, one per frame (set `pdb_folder` in the main block):

    Frames/
        frame0.pdb
        frame500.pdb
        frame1000.pdb
        ...

Each file contains standard fixed-width PDB lines, e.g.:

    ATOM      1  CA1 LIG     1      12.345  10.210   5.678  1.00  0.00
    ATOM      2  CA2 LIG     1      13.101  11.004   6.112  1.00  0.00
    ...

Sample output
-------------
One text file written into the same folder, named with a timestamp:

    Frames/plane_bending_data_2026-09-29_14-05-12.txt

Contents (PDB file name, RMS deviation from plane in Å):

    frame0 0.412837
    frame500 0.398512
    frame1000 0.455901
    ...
"""

import numpy as np
from pathlib import Path
from datetime import datetime
import os

# ---- USER SETTING: folder containing the per-frame PDB files ----
pdb_folder = Path(r"C:\Users\")

def load_pdb(filename):
    """
    Read the x, y, z coordinates of every ATOM/HETATM record in a PDB file.

    Parameters
    ----------
    filename : str or Path
        Path to the PDB file.

    Returns
    -------
    np.ndarray, shape (N, 3)
        Cartesian coordinates (Å) of all N atoms, in file order.
    """
    atoms = []

    with open(filename, "r") as f:
        print(f)  # prints the open file handle (useful for tracking progress)
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # PDB is a fixed-width format; these column slices are standard
                name = line[12:16].strip()   # atom name (not used in the fit)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                #print("Atom {} coords x {}, y {}, z {}".format(name, x,y,z))
                atoms.append((np.array([x, y, z])))

    return np.array(atoms, dtype=float)


def rms_deviation_from_plane(coords):
    """
    RMS distance of a set of points from their least-squares best-fit plane.

    Parameters
    ----------
    coords : np.ndarray, shape (N, 3)
        Atom coordinates. Needs at least 3 points.

    Returns
    -------
    float
        RMS out-of-plane deviation, in the same units as `coords` (Å).
    """
    # Geometric centre of all atoms; the best-fit plane passes through it
    centroid = coords.mean(axis=0)

    # Shift coordinates so the centroid sits at the origin
    centered = coords - centroid

    # Singular value decomposition (equivalent to a principal component
    # analysis of the coordinates). The rows of vh are the principal
    # directions, ordered from most to least spread. The last row, vh[-1],
    # is the direction with the LEAST spread, i.e. the normal of the
    # best-fit plane.
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]

    # Signed perpendicular distance of every atom from the plane
    distances = centered @ normal

    # Root-mean-square of those distances
    rms = np.sqrt(np.mean(distances**2))

    return rms


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
        for name, angle in results:  # `angle` here holds the RMS value
            f.write(f"{name} {angle:.6f} \n")


if __name__ == "__main__":
    os.system("cls")  # clear the console (Windows only)

    results = []

    for pdb_file in pdb_folder.glob("*.pdb"):
        atoms = load_pdb(pdb_file)

        # A plane needs at least 3 points; skip files that don't have them
        if atoms.shape[0] < 3:
            continue

        rms_plane = rms_deviation_from_plane(atoms)
        results.append((pdb_file.stem, rms_plane))

    # One summary file for the whole folder
    write_results(results, pdb_folder, "plane_bending")