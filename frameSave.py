"""
FrameSave.py
============
Author  : Nicholas Bainbridge (nbainbr@student.ubc.ca)
Date    : 02/09/2023
Version : 1.0
Requires: PyMOL (run inside PyMOL, not as a standalone Python script)
 
Description
-----------
PyMOL extension that exports individual frames (states) of a loaded MD
trajectory as separate PDB files. You give a start frame, stop frame and
an interval, and every Nth frame in that range is saved.
 
CONECT records are written for all bonds and the original atom IDs are
kept, so the exported PDBs keep the same atom numbering as the trajectory.
 
Usage
-----
1. Load your topology + trajectory in PyMOL, e.g.
       load system.psf
       load_traj trajectory.dcd, system
2. Drag .py script into pymol
       run FrameSave.py
3. CD into dsave directory
4. Call:
       framesave <filename>, <start>, <stop>, <interval> [, selection]
 
Files are saved to PyMOL's current working directory (check with `pwd`,
change with `cd`).
 
Sample input
------------
Export every 500th frame from frame 0 to 2000:
 
    framesave test, 0, 2000, 500
 
Only save part of the system (e.g. just the peptide, no solvent):
 
    framesave test, 0, 2000, 500, not resn TIP3
 
Sample output
-------------
One PDB per exported frame, named by frame number:
 
    0.pdb, 500.pdb, 1000.pdb, 1500.pdb, 2000.pdb
"""
 
from pymol import cmd, stored
 
 
def framesave(filename, start, stop, interval, selection='(all)', state=0, format='pdb', ref='',ref_state=-1, quiet=1, partial=0):
    """
    Save every `interval`-th state from `start` to `stop` (inclusive) as a PDB.
 
    Parameters
    ----------
    filename : str
        Intended output file prefix (currently unused, see header note).
    start, stop, interval : int or str
        First frame, last frame, and step between saved frames.
        PyMOL passes command-line arguments as strings; they are cast to int.
    selection : str
        PyMOL selection to save. Default '(all)'.
    format : str
        File format passed to cmd.save. Default 'pdb'.
    state, ref, ref_state, quiet, partial :
        Accepted for compatibility with cmd.save's signature; not used.
    """
    # PDB export settings
    cmd.set('pdb_conect_all', 1)    # write CONECT records for all bonds
    cmd.set('pdb_conect_nodup', 0)  # allow duplicate CONECT entries (bond order)
    cmd.set('pdb_retain_ids', 1)    # keep original atom IDs instead of renumbering
 
    # Start at the first requested frame
    x = int(start)
 
    # Iterate over every frame to export, stopping at (and including) the last frame
    while int(x) <= int(stop):
 
        # Output file name: <frame number>.pdb
        newfilename = str(x)+'.pdb'
 
        # Save the selection at state (frame) x
        cmd.save(newfilename, selection, x,format)
 
        # Step to the next frame
        x = x + int(interval)
 
 
# Register `framesave` as a PyMOL command
cmd.extend("framesave", framesave)