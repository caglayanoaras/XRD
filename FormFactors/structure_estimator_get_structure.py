import os, json, ast, random
import pandas                 as     pd
import numpy                  as     np
import matplotlib.pyplot      as     plt
import pymatgen               as     mg

from   XRD_ict                import XRDCalculator
from   pymatgen.core.spectrum import Spectrum
from   pymatgen               import Lattice, Structure
from   diffractogram          import create_perfect_TiAlN_structure, plot_profiles_amplitudes
from   XRD_ict                import XRDCalculator                                    
                                     
plt.rc('font', family='serif')
plt.style.use(['seaborn-paper', 'presentation'])
typ = 'dft'
def error_func(structure, measurement):
    
    calc = XRDCalculator(wavelength="CuKa1", typ = typ)
    XRD = calc.get_pattern(structure,two_theta_range=(30, 105),scaled = True)
    theta_calculated = (np.around(XRD.x, 1))
    two_theta_list = [37.4, 43.5, 63.2, 75.8, 79.8, 95.6 ]
    indexes = []
    for two_theta in two_theta_list:
        indexes.append(np.asscalar(np.argwhere(theta_calculated == two_theta)[0]))

    intensities = XRD.y[indexes]
    return np.square(intensities - measurement).sum()

if __name__ == '__main__':
    TiAlN_perfect      = create_perfect_TiAlN_structure(550, 4.16) #216 atom random structure
    
    D5000 = pd.read_excel('D5000_560.xlsx',sheet_name='DeconvSample')
    measurement = D5000
    
    x = measurement['peak_two_thetas']
    y = measurement['amplitude']
    spec = Spectrum(x,y)
    spec.normalize(mode = 'max', value = 100)
    
    EPOCHS = 400
    UNIT_COST = 0.008
    RISE = 1.025

    
    print('Initial error:{}'.format(error_func(TiAlN_perfect,spec.y)))
    print('-*-*-*-*-*-')

    energies    = []
    structure   = TiAlN_perfect.copy()
    metal_count = len(structure)//2
    length      = len(structure)
    
    for eph in range(EPOCHS):
        energy_current = error_func(structure,spec.y)
        energies.append(energy_current)
        
        predicted_structure = structure.copy()
        # Get a random index
        ind = random.randint(0,length )
        el = predicted_structure[ind].species.formula
        if 'N' in el:
            predicted_structure[ind].species = mg.Composition('V1')
        if 'Ti'in el:
            dice = random.randint(0,1 )
            if dice:
                predicted_structure[ind].species = mg.Composition('V1')
            else:
                predicted_structure[ind].species = mg.Composition('Al1')
        if 'Al'in el:
            dice = random.randint(0,1 )
            if dice:
                predicted_structure[ind].species = mg.Composition('V1')
            else:
                predicted_structure[ind].species = mg.Composition('Ti1')
        if 'Vl'in el:
            if ind < metal_count:
                dice = random.randint(0,1 )
                if dice:
                    predicted_structure[ind].species = mg.Composition('Ti1')
                else:
                    predicted_structure[ind].species = mg.Composition('Al1')
            if ind >= metal_count:
                predicted_structure[ind].species = mg.Composition('N1')
                
        energy_predicted = error_func(predicted_structure,spec.y)

        delta_cost = UNIT_COST*(energy_predicted - energy_current)
        print('Iteration: ', eph)
        print('Energy predicted:{}'.format(energy_predicted))
        print('Energy current:{}'.format(energy_current))
        
        print(np.exp(-delta_cost))
        if np.random.rand() < np.exp(-delta_cost):
            print('switched')
            structure = predicted_structure

        UNIT_COST *= RISE
        print('---')

structure.to(fmt = 'json' , filename = 'best_structure_{}.json'.format(typ))
