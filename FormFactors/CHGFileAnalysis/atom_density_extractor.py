# This script will extract the electron density of a single atom from the huge
# AECCAR file of superlattice. 3D mesh is saved in npy format

from os import listdir
from os.path import isfile, join

import numpy as np
import matplotlib.pyplot as plt

center = np.array([192*2, 192*1, 192*3]) # Center loc of atom
all_layers = np.arange(768)
start      = (center[2] - 96)
stop       = (center[2] +96)
layers     = all_layers[start:stop]

stack_layers = []
for layer in layers:
    edist_line = np.load('layers//layer_{}.npy'.format(layer))
    edist_line_atom = edist_line[(center[0] - 96):(center[0] +96),(center[1] - 96):(center[1] +96)]
    stack_layers.append(edist_line_atom)

    
edist_atom = np.dstack(stack_layers)/(768*768*768)
np.save('edist_atom_{}_{}_{}.npy'.format(center[0],center[1],center[2]), edist_atom)    

plt.imshow(edist_atom[:,:,96])
plt.colorbar()
plt.figure()
plt.imshow(stack_layers[96])
plt.show()

