# This script is for splitting the huge AECCAR files.
# It will create 2D layers from the 3D mesh. 

import numpy as np
import matplotlib.pyplot as plt

mesh = 768

def accar_surf_reader(path, slice_index):
    layer = mesh*mesh
    first_index = layer*slice_index
    final_index = layer*(slice_index+1)
    
    starter = 75
    first_index_in = starter + first_index//5
    final_index_in = starter + final_index//5
    print('First_index_in: ', first_index_in)
    print('Final_index_in: ', final_index_in)    
    big_array = []
    
    with open(path) as f:
        count = 0        
        for line in f:
            count = count +1
            if count >= first_index_in and count <= final_index_in:
                clean_line = line[1:-1].rsplit(" ")
                big_array.extend(clean_line)
            if count > final_index_in:
                break

    initial_deletion = first_index % 5
    final_deletion   = 5 - final_index %5
    print('initial_deletion: ', initial_deletion)
    print('final_deletion: ', final_deletion)
    big_array_float = [float(i) for i in big_array]
    surface = np.array(big_array_float)[initial_deletion:-final_deletion]
    
    
    return surface.reshape(mesh,mesh)

for layer in range(740,767):
    path0 = 'AECCAR0'
    surf0 = accar_surf_reader(path0,layer)

    path2 = 'AECCAR2'
    surf2 = accar_surf_reader(path2,layer)
    np.save('layer_{}.npy'.format(layer), surf0+surf2)
