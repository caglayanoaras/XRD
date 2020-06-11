import os, json, ast, random
import pandas                 as     pd
import numpy                  as     np
import matplotlib.pyplot      as     plt

from   XRD_ict                import XRDCalculator
from   pymatgen.core.spectrum import Spectrum
from   pymatgen               import Lattice, Structure, Molecule

plt.rc('font', family='serif')
plt.style.use(['seaborn-paper', 'presentation'])

def create_perfect_TiAlN_structure(number_of_atoms, a, shuffle = True):
    ''' Attributes:
        number_of_lattice > maximum Total number of atoms in structure
        a                 > lattice parameter
        Returns:
        Pymatgen structure TiAlN. Creates possible largest square supercell '''

    number_of_lattice = number_of_atoms //8
    
    lattice = Lattice.cubic(a = a)
    
    matrix_coef = int(number_of_lattice**(1/3))
    
    coords = [[0   , 0   , 0   ],
              [0.5 , 0.5 , 0   ],
              [0.5 , 0   , 0.5 ],
              [0   , 0.5 , 0.5 ],
              [0   , 0   , 0.5 ],
              [0.5 , 0   , 0   ],
              [0   , 0.5 , 0   ],
              [0.5 , 0.5 , 0.5 ],
              ]
    matrix = [[matrix_coef,0,0],
              [0,matrix_coef,0],
              [0,0,matrix_coef]]
    
    struct = Structure(lattice, ["Ti", "Ti", "Al", "Al","N", "N", "N", "N"], coords)
    struct.make_supercell(matrix)
    
    if shuffle == True:
        
        site_size = len(struct)
        element_number = int(site_size/4)
        Al_Ti_array = ['Al'] * element_number + ['Ti'] * element_number
        random.shuffle(Al_Ti_array)
        for e,i in enumerate(Al_Ti_array):
          struct.replace(e,i)
                
    return struct

def plot_profiles_amplitudes(*args, structure = None,labellist = None):
    ''' Attributes:
        args      > list of different measurements in pandas form
        structure > pymatgen structure
        Returns:
        plot of calculated versus measurement patterns
    '''
    calc = XRDCalculator(wavelength="CuKa1", typ = 'dft')
    calc2= XRDCalculator(wavelength="CuKa1", typ = 'ict')
    XRD = calc.get_pattern(structure,two_theta_range=(15, 105),scaled = True)
    XRD2 = calc2.get_pattern(structure,two_theta_range=(15, 105),scaled = True)
    
##    for i in range(0,len(XRD.hkls)):
##        print(XRD.hkls[i])
##        print(XRD.y[i])    
    fig, ax = plt.subplots()
    ax.set_xlim([9, 105])
    ax.set_ylim([0, 130])
    
    ax.scatter(XRD.x, np.multiply(XRD.y, 1.0), marker='X',color = 'b',label = 'Ab initio form factors')
    ax.scatter(XRD2.x, np.multiply(XRD2.y, 1.0), marker='X',color = 'r',label = 'Ict form factors')
    for i,j in zip(XRD.x,np.multiply(XRD.y, 1.0)):
        ax.plot([i,i],[0,j],linestyle = '--', color = 'b' ,alpha = 0.5)
    for i,j in zip(XRD2.x,np.multiply(XRD2.y, 1.0)):
        ax.plot([i,i],[0,j],linestyle = '--', color = 'b' ,alpha = 0.5)        
    for c,i in enumerate(args):
        if c == 0:
            color = 'k'
        if c == 1:
            color = 'r'
        if c == 2:
            color = 'b'
        x = i['peak_two_thetas']
        y = i['amplitude']
        spec = Spectrum(x,y)
        spec.normalize(mode = 'max', value = 100)
        ax.scatter(spec.x, spec.y, marker='X',color = color,label = labellist[c])
        for i,j in zip(spec.x,spec.y):
            ax.plot([i,i],[0,j],linestyle = '--', color = color ,alpha = 0.5)

    if labellist != None:
        plt.legend()
        plt.grid()
        
    plt.xlabel('2\u03B8')
    plt.ylabel('Intensity (a.u.)')
    plt.show()
    
if __name__ == '__main__':
    
    D5000 = pd.read_excel('GADDS_asdep.xlsx',sheet_name='DeconvSample')
    TiAlN_perfect      = create_perfect_TiAlN_structure(700, 4.16)
       
    plot_profiles_amplitudes(D5000, structure = TiAlN_perfect,labellist = ['Measurement'])
