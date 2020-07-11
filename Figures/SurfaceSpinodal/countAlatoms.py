import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from   pymatgen          import Lattice, Structure, Molecule
import os, json, ast, random
from molecularDynamics import create_lammps_structure_file

# Random 
alcont_y = np.array([3,3,3,2,2,3,2,2,2,3]) /5 *100
alcont_z = np.array([3,2,1,3,2,3,3,2,4,2]) /5 *100

alcont_11      = np.array([5,4,7,4,5,5,4,7,4]) /10 *100
alcont_minus11 = np.array([6,6,5,5,3,6,6,5,5,3]) /10 *100

alcount = np.append(alcont_y, alcont_z)
alcount_diag = np.append(alcont_11, alcont_minus11)


# spinodal 
alcont_y_s = np.array([1,1,3,4,5,5,2,1,2,2]) /5 *100
alcont_z_s = np.array([2,2,1,1,5,5,1,2,4,2]) /5 *100

alcont_11_s      = np.array([6,5,6,5,3,6,5,6,5]) /10 *100
alcont_minus11_s = np.array([6,6,5,5,3,6,6,5,5,3]) /10 *100

alcount_s = np.append(alcont_y_s, alcont_z_s)
alcount_diag_s = np.append(alcont_11_s, alcont_minus11_s)

print(np.std(alcount))
print(np.std(alcount_diag))
print('*-*-*-*-')
print(np.std(alcount_s))
print(np.std(alcount_diag_s))
