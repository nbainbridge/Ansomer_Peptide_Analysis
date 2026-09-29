"""
DihedralAngle.py
================
Author  : Nicholas Bainbridge (nbainbr@student.ubc.ca)
Requires: Python 3, NumPy
 
Description
-----------
Batch-calculates dihedral (torsion) angles, e.g. backbone phi/psi, from a
set of MD snapshot PDB files.
 
Each dihedral is defined by four atom names (A-B-C-D) and is the angle
between the A-B-C plane and the B-C-D plane, looking down the B-C bond.
Values are signed and returned in degrees in the range -180 to +180.
 
The script walks every subfolder of `base_dir` (e.g. one subfolder per
simulation/condition), calculates each dihedral in `dihedral_sets` for
every PDB in that subfolder, and writes one results file per dihedral
into that subfolder.
 
Atoms are looked up by their PDB atom name (columns 13-16), so atom names
must be unique within each file.
 
Sample input
------------
    Base_dir/                     <- base_dir
        Water/                   <- one subfolder per condition
            frame0.pdb
            frame500.pdb
            ...
        DMSO/
            frame0.pdb
            ...
 
Dihedrals to measure (set in the main block):
 
    dihedral_sets = [
        ("C7", "N6", "C8", "C14", "Dap_phi"),
        ("N6", "C8", "C14", "N8", "Dap_psi"),
    ]
 
Sample output
-------------
Files written into each subfolder, named <label>_data_<timestamp>.txt:
 
    Base_dir/Water/phi_data_2026-09-29_14-05-12.txt
    Base_dir/Water/psi_data_2026-09-29_14-05-12.txt
 
Contents (PDB file name, dihedral angle in degrees):
 
    frame0 -72.418305
    frame500 -65.902114
    frame1000 -151.337820
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
base_dir = Path(r"E:\")

# ---- USER SETTING: dihedrals to measure ----
# Each entry: (atom1, atom2, atom3, atom4, label used in the output file name)
dihedral_sets = [
    ("C7", "N6", "C8", "C14", "Dap_phi"),
    ("N6", "C8", "C14", "N8", "Dap_psi"),
]


def write_results(results, output_folder, label):
    """
    Write (name, value) pairs to a timestamped text file.
 
    Output file: <output_folder>/<label>_data_<YYYY-mm-dd_HH-MM-SS>.txt
    One line per result: "<name> <value to 6 d.p.>"
    """
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
    output_file = Path(output_folder) / f"{label}_data_{timestamp}.txt"
    with open(output_file, "w") as f:
        for name, distance in results:  # `distance` here holds the angle
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
 
 
def dihedral(atoms, atom1, atom2, atom3, atom4):
    """
    Signed dihedral angle atom1-atom2-atom3-atom4, in degrees (-180 to 180).
 
    Parameters
    ----------
    atoms : dict
        {atom_name: coordinates} as returned by load_pdb().
    atom1, atom2, atom3, atom4 : str
        Atom names defining the dihedral. The rotation is about atom2-atom3.
 
    Raises
    ------
    ValueError
        If any of the four atoms is not present in `atoms`.
    """
    # Make sure all four atoms exist in this frame
    if atom1 not in atoms:
        raise ValueError(f"{atom1} not found")
 
    if atom2 not in atoms:
        raise ValueError(f"{atom2} not found")
 
    if atom3 not in atoms:
        raise ValueError(f"{atom3} not found")
 
    if atom4 not in atoms:
        raise ValueError(f"{atom4} not found")
 
    p0 = atoms[atom1]
    p1 = atoms[atom2]
    p2 = atoms[atom3]
    p3 = atoms[atom4]
 
    # Bond vectors along the chain
    b0 = -(p1 - p0)   # atom2 -> atom1
    b1 = p2 - p1      # atom2 -> atom3 (the central bond / rotation axis)
    b2 = p3 - p2      # atom3 -> atom4
 
    # Normalize b1 so the projections below are simple dot products
    b1 /= np.linalg.norm(b1)
 
    # Components of b0 and b2 perpendicular to b1 (i.e. projected onto the
    # plane perpendicular to the central bond)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
 
    # Angle between v and w in that plane. arctan2 gives the correct sign
    # and quadrant (-180 to 180) without numerical issues from arccos.
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
 
    return np.degrees(np.arctan2(y, x))
 
 
if __name__ == "__main__":
    os.system("cls")  # clear the console (Windows only)

 
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
 
        # One results file per dihedral, per subfolder
        for atom1, atom2, atom3, atom4, label in dihedral_sets:
 
            results = []
            for pdb_file in pdb_files:
                atoms = load_pdb(pdb_file)
                try:
                    angle = dihedral(atoms,atom1,atom2,atom3,atom4)
                    results.append((pdb_file.stem, angle))
                except ValueError as e:
                    # Report frames missing an atom and carry on
                    print(f"{pdb_file.name}: {e}")
 
            if results:
                write_results(results, subdir, label)
 
            print(f"Processed {subdir.name}")