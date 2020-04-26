#miscellaneous

from os import listdir
from os.path import isfile, join

import numpy as np
import matplotlib.pyplot as plt

def show_layer(layer = None, log = True):
    
    edist = np.load('layers//layer_{}.npy'.format(layer))
    
    if log ==True:
        plt.imshow(np.log(edist))
    else:
        plt.imshow(edist)
        
    plt.xticks([])
    plt.yticks([])
    plt.colorbar()
    plt.show()

def checksum(layers = np.arange(0,768)):
    #check the sum of all layers
    sumlist = np.zeros(layers.size)
    
    for c,layer in enumerate(layers):
        
        edist = np.load('layers//layer_{}.npy'.format(layer))

        sumlist[c] = edist.sum()

    return sumlist
    


if __name__ == '__main__':

    show_layer(layer = 0, log = True)

    sumlist = checksum()
    electron = sumlist.sum() / (768*768*768)
    print('number of electrons should be 784: ', electron)
