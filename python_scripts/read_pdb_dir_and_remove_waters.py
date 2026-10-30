# AD 10.10.2026
import argparse
from requests import get
import os
import glob
import MDAnalysis as md
####
"""
here we read all pdb files in a directory and 
rewrite them as their plain 4-letetr PDB name
and also remove any waters.
"""
parser = argparse.ArgumentParser()
parser.add_argument("pdbdir",type=str,help="directory containing PDB files with names aligned_PDBID.pdb")
args = parser.parse_args()
####
print("# USING DIR",args.pdbdir)
files = glob.glob(args.pdbdir+"/*.pdb")
print("# CONTAINS",len(files),"FILES")
print("# FILE 0:",files[0])

ind = int()
for line in files:
    try:
        u = md.Universe(line)
        prot = u.select_atoms("protein and not resname HOH")
        prot.atoms.write(f"fixed_{line.split('_')[-1]}")
        ind += 1
        # get rid of old pdb
    except Exception as exc:
        print("#FAILED DURING:",line,"AT INDEX:",ind)
        print("#REASON:",exc)
        continue
########################## COMPLETE
print("# CLEAN COMPLETE: KYRIRE ELESON")

