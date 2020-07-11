import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from   pymatgen          import Lattice, Structure, Molecule
import os, json, ast, random
from molecularDynamics import create_lammps_structure_file
import sys
sys.exit()
def create_TiAlN_structure( shuffle = True):
    ''' Attributes:
        number_of_lattice > maximum Total number of atoms in structure
        a                 > lattice parameter
        Returns:
        Pymatgen structure TiAlN. Creates possible largest square supercell '''
    lattice = Lattice.cubic(a = 4.16)
    coords = [[0   , 0   , 0   ],
              [0.5 , 0.5 , 0   ],
              [0.5 , 0   , 0.5 ],
              [0   , 0.5 , 0.5 ],
              [0   , 0   , 0.5 ],
              [0.5 , 0   , 0   ],
              [0   , 0.5 , 0   ],
              [0.5 , 0.5 , 0.5 ],
              ]
    matrix = [[1,0,0],
              [0,5,0],
              [0,0,5]]
    
    struct = Structure(lattice, ["Ti", "Ti", "Al", "Al","N", "N", "N", "N"], coords)
    struct.make_supercell(matrix)
    site_size = len(struct)
    element_number = int(site_size/4)
    Al_Ti_array = ['Al'] * element_number + ['Ti'] * element_number
    
    if shuffle == True:
        random.shuffle(Al_Ti_array)
        for e,i in enumerate(Al_Ti_array):
          struct.replace(e,i)
    return struct

structure = create_TiAlN_structure( shuffle = True)

metal_indexes_of_front =[]
for e,site in enumerate(structure):
    if site.x == 0 and (site.species.formula == 'Al1' or site.species.formula == 'Ti1'):
        metal_indexes_of_front.append(e)
metal_number = len(metal_indexes_of_front)

Al_Ti_array = ['Al', 'Ti'] * int(metal_number/2) 
random.shuffle(Al_Ti_array)

for e,i in enumerate(metal_indexes_of_front):
    structure[i].species = Al_Ti_array[e]

    
structure.to(filename = 'SurfaceRandom.json')

ALN_list = [10,11,12,13,14,35,36,37,38,39] + [2,7,12,17,22,27,32,37,42,47]
anti_ALN_list = [i for i in range(0,len(metal_indexes_of_front)) if i not in ALN_list]

Al_Ti_array_for_remainig = ['Al'] * 7 + ['Ti']* 25
random.shuffle(Al_Ti_array_for_remainig)

structure_spi = structure.copy()
for i in ALN_list:
    structure_spi[metal_indexes_of_front[i]].species = 'Al'
for e,i in enumerate(anti_ALN_list):
    structure_spi[metal_indexes_of_front[i]].species = Al_Ti_array_for_remainig[e]
structure_spi.to(filename = 'SurfaceSpinodal.json')
