"""
twoAtomDistance.py
==================
Author  : Nicholas Bainbridge (nbainbr@student.ubc.ca)
Requires: Python 3, NumPy
 
Description
-----------
Batch-calculates interatomic distances (e.g. NH-NH distances or H-bond
donor/acceptor distances) for a list of atom pairs across a set of MD
snapshot PDB files.
 
The distance is the straight-line (Euclidean) distance between the two
atoms, in Å.
 
The script walks every subfolder of `base_dir` (e.g. one subfolder per
simulation/condition), measures every pair in `atom_pairs` for every PDB
in that subfolder, and writes one results file per atom pair into that
subfolder. Frames missing one of the atoms are skipped and reported.
 
Atoms are looked up by their PDB atom name (columns 13-16), so atom names
must be unique within each file.
 
Sample input
------------
    LowE_Aq/                     <- base_dir
        Water/                   <- one subfolder per condition
            frame0.pdb
            frame500.pdb
            ...
        DMSO/
            frame0.pdb
            ...
 
Atom pairs to measure (set in the main block):
 
    atom_pairs = [
        ("H1", "H2"),
        ("H2", "O7"),
        ...
    ]
 
Sample output
-------------
One file per atom pair written into each subfolder, named
<atom1>_<atom2>_data_<timestamp>.txt:
 
    Base_dir/Water/H1_H2_data_2026-09-29_14-05-12.txt
    Base_dir/Water/H2_O7_data_2026-09-29_14-05-12.txt
    ...
 
Contents (PDB file name, distance in Å):
 
    frame0 2.714392
    frame500 3.108825
    frame1000 2.455017
    ...
 
Console:
 
    Processed Water
    ...
"""
 
import numpy as np
from pathlib import Path
from datetime import datetime
import os
 
 
def write_results(results, output_folder, atom1, atom2):
    """
    Write (name, distance) pairs to a timestamped text file.
 
    Output file: <output_folder>/<atom1>_<atom2>_data_<YYYY-mm-dd_HH-MM-SS>.txt
    One line per result: "<name> <distance to 6 d.p.>"
    """
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
    output_file = Path(output_folder) / f"{atom1}_{atom2}_data_{timestamp}.txt"
    with open(output_file, "w") as f:
        for name, distance in results:
            f.write(f"{name} {distance:.6f} \n")
 
 
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
                atom_name = line[12:16].strip()
 
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
 
                atoms[atom_name] = np.array([x, y, z])
 
    return atoms
 
 
def twoPointDist(atoms, atom1, atom2):
    """
    Euclidean distance between two named atoms, in the PDB's units (Å).
 
    Raises
    ------
    ValueError
        If either atom is not present in `atoms`.
    """
    if atom1 not in atoms:
        raise ValueError(f"{atom1} not found")
 
    if atom2 not in atoms:
        raise ValueError(f"{atom2} not found")
 
    # |r1 - r2|
    return np.linalg.norm(atoms[atom1] - atoms[atom2])
 
 
if __name__ == "__main__":
    os.system("cls")  # clear the console (Windows only)
 
    # ---- USER SETTINGS: data folder and atom pairs to measure ----
    # Swap between the blocks below by commenting/uncommenting.
 
    # LowEAq/DMSO NH-NH Distances (1-7) and Internal H-Bonding Network (8-10)
    base_dir = Path(r"E:\Desktop\LowE_Aq")
    atom_pairs = [
        ("H1", "H2"),
        ("H2", "H3"),
        ("H2", "H4"),
        ("H3", "H4"),
        ("H14", "H5"),
        ("H14", "H6"),
        ("H5", "H6"),
        ("H2", "O7"),
        ("H4", "O7"),
        ("H14", "O4"),
    ]
 
    # ND1115 Atom Distances
    # base_dir = Path(r"E:\Desktop\ND11_15")
    # atom_pairs = [
    #     ("H12", "H13"),
    #     ("H9", "H12"),
    #     ("H9", "H13"),
    #     ("H1", "H13"),
    #     ("H2", "H3"),
    #     ("H3", "H8"),
    #     ("H7", "H8"),
    # ]
 
    # Loop over each subfolder in base_dir
    for subdir in base_dir.iterdir():
        if not subdir.is_dir():
            continue  # skip files in base directory
        pdb_files = list(subdir.glob("*.pdb"))
 
 
        # Find PDB files inside this subfolder
        pdb_files = list(subdir.glob("*.pdb"))
        if not pdb_files:
            continue  # skip folders with no pdbs
 
        if not pdb_files:
            continue
 
        # One results file per atom pair, per subfolder
        for atom1, atom2 in atom_pairs:
 
            results = []
            for pdb_file in pdb_files:
                atoms = load_pdb(pdb_file)
                try:
                    distance = twoPointDist(atoms, atom1, atom2)
                    results.append((pdb_file.stem, distance))
                except ValueError as e:
                    # Report frames missing an atom and carry on
                    print(f"{pdb_file.name}: {e}")
 
            if results:
                write_results(results, subdir, atom1, atom2)
 
            print(f"Processed {subdir.name}")
 