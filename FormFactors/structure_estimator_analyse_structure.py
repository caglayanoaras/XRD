import pymatgen          as mg
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
from   XRD_ict           import XRDCalculator
from   diffractogram     import plot_profiles_amplitudes

best_structure = mg.Structure.from_file('best_structure_dft.json')
D5000 = pd.read_excel('D5000_560.xlsx',sheet_name='DeconvSample')

specie_dict = {'element':best_structure.species}
specie_df = pd.DataFrame(specie_dict)

print(specie_df['element'].value_counts())
plot_profiles_amplitudes(D5000, structure = best_structure,labellist = ['Experimental'])
