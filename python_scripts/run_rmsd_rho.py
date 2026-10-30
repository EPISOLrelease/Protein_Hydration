# AD 10.11.2026 ICXC NIKA
"""
read in three pdbfiles:
1) 3DRISM placed waters
2) SuperWater placed waters
3) true waters

align both structures to the true
and then write out and take RMSD
and PREC and RECALL

select 3DRISM above a certain RHO threshold
from the tempfactors. 
sort the resulting positions based on density
then remove any waters that are within RADII
with higest dens first.
DO NOT update the positions, just remove the overlaps
"""
import numpy as np 
from scipy import ndimage
## personal packages
from download_rcsb import *
from epipy import *
from utils import *
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
###
import warnings

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
#import pymol 
import time
import numpy as np
import argparse
import MDAnalysis as md
from MDAnalysis.lib.distances import capped_distance
import MDAnalysis.transformations as mdt
import getpass #DOES NOT WORK ON CHTC
import datetime
import os
# ----------------------------------------------------
#                   CLUSTER ALG
# ----------------------------------------------------

def cluster_and_return_coords(coord_list,radii,use_mean):
    """ 
    read in a series of coordinates 
    and cluster the coordinates in descending order
    return a list of cluster centers grouped by mean

    >WARNING:
    assumes coords are already sorted
    """
    assigned = np.zeros(len(coord_list),dtype=bool)
    out_list = []
    for i in range(len(coord_list)):
        #print('shape',i,assigned.shape)
        #print('cord shape',i,coord_list.shape)
        if assigned[i]:
            continue
        else:
            distances = np.linalg.norm(coord_list[i] - coord_list,axis=1)
            mask = (distances < radii) & (~assigned)
            #
            cluster = coord_list[mask]
            assigned[mask] = True
            #
            #coord_list = coord_list[~mask]
            #print("reduced to length",len(coord_list))
            if use_mean:
                out_list.append(cluster.mean(axis=0))
            else:
                out_list.append(coord_list[i])

    return np.array(out_list)
# ----------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("name",help="AF PDBID",type=str)
parser.add_argument("DENS",help="density gt")#,type=str)
args = parser.parse_args()
print(f"running on: {getpass.getuser()}")
print(f"# PARSED: {args.name}")
#################
# ---------------------------------------------------- PARAMS
GT = float(args.DENS) # place at regions with density greater than this value
N_STEPS = 10
TEMP = 300
USE_NO_MASK = True
SIGMA = 0.5 # uncertainty distance in 3drism grid
RADII = 2.2#1.52
USE_MEAN = False # use mean of cluster centroid for new position 
# -------------- END PARAMS
print("# USING NO_MASK",USE_NO_MASK)
print("# USING TEMP",GT)
print("# USING N_STEPS",GT)
print("# USING RHO GT",GT)
# ------------------------------------------- rewrite superwater to match old
md.Universe(f"aligned_combined_{args.name}.pdb").atoms.write(f"combined_{args.name}.pdb")

# ------------------------------------------------- CLUSTER 3DRISM

total_placed  = md.Universe(f'aligned_{args.name}.pdb')#.select_atoms('resname HOH').positions
#
rism_positions = total_placed.select_atoms("resname HOH")
rism_densities = rism_positions.atoms.tempfactors
# NEW 
rism_inds = np.where(rism_densities > GT )
rism_positions = total_placed.select_atoms("resname HOH").positions

rism_densities = rism_densities[rism_inds]
rism_positions = rism_positions[rism_inds]


# order positions based on the highest densities
# numpy argsort deafults to ascending ... lol
rism_sorted = rism_positions[np.argsort(rism_densities)[::-1]]
print("# STARTING N_POSITIONS:",len(rism_sorted))

new_positions = cluster_and_return_coords(coord_list=rism_sorted,radii=RADII,use_mean=USE_MEAN)
#---------------------------------------------------- CREATE NEW PROT+WATER
rism_prot = total_placed.select_atoms("protein")

n_atoms = len(new_positions)
seg_index = np.zeros(n_atoms)
res_index = seg_index+1
wot = md.Universe.empty(n_atoms=n_atoms,n_residues=n_atoms,
                        atom_resindex=res_index,residue_segindex=seg_index,trajectory=True)

wot.add_TopologyAttr('name',['O']*n_atoms)
wot.add_TopologyAttr('resname',['HOH']*n_atoms)
wot.add_TopologyAttr('segid',['HOH'])
wot.add_TopologyAttr('resid',[i for i in range(n_atoms)])
wot.add_TopologyAttr('ids',[i for i in range(n_atoms)])
# write out
wot.atoms.positions = new_positions
combined = md.Merge(rism_prot.atoms,wot.atoms)
#combined.atoms.write("{args.name}_rism_placed_waters.pdb")
combined.dimensions = rism_prot.dimensions
#combined.add_TopologyAttr('tempfactors',np.append(np.zeros(len(prot.atoms)),dens_out))

combined.atoms.write(f"{args.name}_rism_placed_waters.pdb")

print("# CLUSTERED N_POSITIONS:",len(new_positions))

#-------------------------------------------------- ALIGN
# placed from 3DRISM
#placed  = md.Universe('{args.name}_rism_placed_waters.pdb')#.select_atoms('resname HOH').positions
placed  = md.Universe(f"{args.name}_rism_placed_waters.pdb")

# (f'aligned_{args.name}.pdb')
#.select_atoms('resname HOH').positions
# true waters
tru = md.Universe(f'full_{args.name}.pdb')
# waters placed by superwater
supr = md.Universe(f"aligned_combined_{args.name}.pdb")
#
num_orig_wat = len(tru.select_atoms("resname HOH"))
num_placed = len(placed.select_atoms("resname HOH"))
num_supr = len(supr.select_atoms("resname HOH"))

## set box sizes to be the same as precautionary

tru.dimensions = placed.dimensions
## select the backbone CA positions
# ----------------------------------------------------------------- START ALIGN
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored

pymol.finish_launching()

# load into pymol for alignment
cmd.load(f"full_{args.name}.pdb")
cmd.load(f"{args.name}_rism_placed_waters.pdb")
cmd.load(f"combined_{args.name}.pdb")
print("# LOADED INTO PYMOL")
# Now align the RISM, and the supoerwater placed 
# WARNING! WILL FAIL FOR NUCLEOTIDE. but I am only doing prot. so this is a failsafe
cmd.align(f"polymer and name CA and ({args.name}_rism_placed_waters)",f"polymer and name CA and (full_{args.name})",quiet=0,object="aln",reset=1)
cmd.align(f"polymer and name CA and (combined_{args.name})",f"polymer and name CA and (full_{args.name})",quiet=0,object="aln",reset=1)
print("# ALIGNED, SAVING...")
cmd.save(f"aligned_{args.name}.pdb",f"{args.name}_rism_placed_waters")
cmd.save(f"aligned_combined_{args.name}.pdb",f"combined_{args.name}")
print("# SAVED")
cmd.reinitialize()
#-------------------------------------------------- END ALIGN


# ------------------------------------------------- SHELLS

tru = md.Universe(f'full_{args.name}.pdb')
# 3DRISM aligned
placed = md.Universe(f"aligned_{args.name}.pdb")
# SuperWater aligned
suprwtr = md.Universe(f"aligned_combined_{args.name}.pdb")

# for each solvation shell, around each residue count number of waters
print("# RESNAME LOW HIGH MU-NUMBER-OF-WATERS STDEV-NUMBER-OF-WATERS VAR-NUMBER-OF-WATERS")
SHELLS = [(0.0,1.0),(1.0,2.0),(2.0,3.0),(3.0,4.0)]
#
for u_name,u_i in zip(["TRU","3DP","SWP"],[tru,placed,suprwtr]):
    print("#",17*"-",u_name,18*"-")
    for resname_i in np.unique(u_i.select_atoms("protein").resnames):
        for (low_i,hi_i) in SHELLS:
            to_go = []
            for name_i in u_i.select_atoms(f"resname {resname_i}").resids:
                N_i = len(u_i.select_atoms(f"resname HOH and (sphlayer {low_i} {hi_i} resid {name_i})"))
                to_go.append(N_i)
            to_go = np.array(to_go)
            print(f"{u_name} {resname_i} {low_i} {hi_i} {to_go.mean():0.3f} {to_go.std():0.3f} {to_go.var():0.3f}")
# ------------------------------------------------- END SHELLS


# ------------------------------------------------- SHELLS AROUND PLACED WATER
print("# SHELL AROUND PLACED WATERS")
print("# NAME SHELLNAME LOW HI COUNT")
SHELLS = [(0.0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.5,3.0),(3.5,4.0),(4.5,5.0),(5.0,5.5),(5.5,6.0)]
#
#
tru_wat_pos = tru.select_atoms("resname HOH").positions # select all waters
# SELECT WATER IN SPHERICAL SHELLS AROUND TRU WAT
for u_name,u_i in zip(["3DP","SWP"],[placed,suprwtr]):
    # make a new dictionary
    td = dict(zip([f'shell_{i}' for i in range(len(SHELLS))],[[i,int(),int()] for i in SHELLS]))
    print("#",17*"-",u_name,18*"-")
    if u_name == '3DP':
        sigma = SIGMA #1.52 # smearing function
    else:
        sigma = 0.0
    for i in u_i.select_atoms("resname HOH"):
        for tru_wat_pos_i in tru_wat_pos:
            per_wat_id = int() # for each water. if there are too many waters
            for key,((lo,hi),c,c1) in td.items():
                # take distance
                x = np.linalg.norm(i.position - tru_wat_pos_i) - sigma
                tmp_dist = np.piecewise(x,[x<=0,x>0],[0.1,x])
                if lo < tmp_dist <= hi:
                    td[key][1] += 1
                    per_wat_id += 1 # keep track of surrounding wat

                td[key][2] += per_wat_id
    # --------------- write to stdout 
    for key,((lo,hi),c,c1) in td.items():
        print(f"{u_name} {key} {lo} {hi} {c:0.3f} {c1:0.3f}")
# ------------------------------------------------- END SHELLS AROUND PLACED WATER

# ------------------------------------------------- PREC and RECALL
print("#\n# SHELL AROUND CRYSTAL WATERS")
print("# NAME SHELLNAME LOW HI COUNT")
SHELLS = [(0.0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.5,3.0),(3.5,4.0),(4.5,5.0),(5.0,5.5),(5.5,6.0)]
#
#
tru = md.Universe(f'full_{args.name}.pdb')
tru_wat_pos = tru.select_atoms("resname HOH").positions # select all waters
# SELECT WATER IN SPHERICAL SHELLS AROUND TRU WAT
for u_name,u_i in zip(["3DP","SWP"],[placed,suprwtr]):
    # make a new dictionary
    print("#",17*"-",u_name,18*"-")
    td = dict(zip([f'shell_{i}' for i in range(len(SHELLS))],[[i,int(),int()] for i in SHELLS]))
    if u_name == '3DP':
        sigma = SIGMA # subtract uncertainty from 3drism grid
    else:
        sigma = 0.0
    for tru_wat_pos_i in tru_wat_pos:
        # now select the placed waters
        for i in u_i.select_atoms("resname HOH"):
            per_wat_id = int() # for each water.
            for key,((lo,hi),c,c1) in td.items():
                # take distance
                x = np.linalg.norm(i.position - tru_wat_pos_i) - sigma
                tmp_dist = np.piecewise(x,[x<=0,x>0],[0.1,x])
                if lo < tmp_dist <= hi:
                    td[key][1] += 1
                    per_wat_id += 1 # keep track of surrounding wat

                td[key][2] += per_wat_id
    # --------------- write to stdout 
    for key,((lo,hi),c,c1) in td.items():
        print(f"CWAT {u_name} {key} {lo} {hi} {c:0.3f} {c1:0.3f}")

# ------------------------------------------------- END PREC and RECALL

# ------------------------------------------------- RMSD 
# select the protein atoms from the true PDB
tru_prot = tru.select_atoms("protein") # make atomgroup
print("# ABOUT TO COMPUTE DISTANCES")
#------------------- PASS NAME TO WATERS ONLY 
placed  = placed.select_atoms("resname HOH")
suprwtr = suprwtr.select_atoms("resname HOH")

# merge two selections together and select the waters
placed = md.Merge(tru_prot.atoms,placed.atoms).select_atoms("resname HOH").positions
suprwtr = md.Merge(tru_prot.atoms,suprwtr.atoms).select_atoms("resname HOH").positions

# select two copies
supr_tru = tru.select_atoms("resname HOH").positions
rism_tru = tru.select_atoms("resname HOH").positions

def match_points_and_get_distances(data1, data2):
    """
    from superwater .git
    """
    dist_matrix = cdist(data1, data2)
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    min_distances = np.array([dist_matrix[i, j] for i, j in zip(row_ind, col_ind)])
    return min_distances

print("# ================== START 3DRISM")
tt = match_points_and_get_distances(placed,rism_tru) 

print("ORIG",num_orig_wat)
print("PLACED_R",num_placed)

print(f"LSA RISM {np.mean(tt):0.3f}")
print(f"MUSQ_LSA_RISM {np.sqrt(np.mean(tt**2)):0.3f}")
x = match_points_and_get_distances(placed,rism_tru) - SIGMA
tt = np.where(x <= 0 , 0, x)
print(f"MUSQ_LSA_RISM_W_SIGMA {np.sqrt(np.mean(tt**2)):0.3f}")
tt = cdist(placed,rism_tru)
print(f"R_NEAREST_1 {tt.min(axis=0).mean():0.3f}")
print(f"R_NEAREST_2 {tt.min(axis=1).mean():0.3f}")

tt = match_points_and_get_distances(placed,placed) 
print(f"R_TRU_TRU {np.mean(tt):0.3f}")

# ======================= PREC and RECALL no double count
tt = cdist(placed,rism_tru)
print(f"R_FP_05 {np.float64(tt.min(axis=1) <= 0.5+SIGMA).mean() :0.3f}")
print(f"R_FP_1  {np.float64(tt.min(axis=1) <= 1+SIGMA).mean() :0.3f}")

print(f"R_TP_05 {np.float64(tt.min(axis=0) <= 0.5+SIGMA).mean() :0.3f}")
print(f"R_TP_1  {np.float64(tt.min(axis=0) <= 1+SIGMA).mean() :0.3f}")
# 
print("==================")
tt = cdist(placed,placed)
tt = np.where(tt == 0, 99 ,tt)
print(f"R_PLACED_PLACED {tt.min(axis=0).mean():0.3f}")
print(f"R_PLACED_PLACED_2 {tt.min(axis=1).mean():0.3f}")
print("================== END 3DRISM")

print("================== START SUPERWATER")
print("ORIG",num_orig_wat)
print("PLACED_S",num_supr)
tt = match_points_and_get_distances(suprwtr,supr_tru)
print(f"LSA SUPRWTR {np.mean(tt):0.3f}")
print(f"MUSQ_LSA_SUPRWTR {np.sqrt(np.mean(tt**2)):0.3f}")
x = match_points_and_get_distances(suprwtr,supr_tru) - SIGMA
tt = np.where(x <= 0 , 0, x)
print(f"MUSQ_LSA_SUPRWTR_W_SIGMA {np.sqrt(np.mean(tt**2)):0.3f}")

print("==================")
tt = cdist(suprwtr,supr_tru)
print(f"S_NEAREST_1 {tt.min(axis=0).mean():0.3f}")
print(f"S_NEAREST_2 {tt.min(axis=1).mean():0.3f}")
tt = match_points_and_get_distances(suprwtr,supr_tru)
print(f"S_TRU_TRU {np.mean(tt):0.3f}")

# ======================= PREC and RECALL no double count
tt = cdist(suprwtr,supr_tru)
print(f"S_FP_05 {np.float64(tt.min(axis=1) <= 0.5).mean() :0.3f}")
print(f"S_FP_1  {np.float64(tt.min(axis=1) <= 1).mean() :0.3f}")

print(f"S_TP_05 {np.float64(tt.min(axis=0) <= 0.5).mean() :0.3f}")
print(f"S_TP_1  {np.float64(tt.min(axis=0) <= 1).mean()  :0.3f}")
print("==================")
tt = cdist(suprwtr,suprwtr)
tt = np.where(tt == 0,99,tt)
print(f"S_PLACED_PLACED_1 {tt.min(axis=0).min():0.3f}")
print(f"S_PLACED_PLACED_2 {tt.min(axis=1).min():0.3f}")
print("================== END SUPERWATER")
###################################
print('# Run completed, deleting solute, topology and dump files')
print("# hope you have a good day!")

