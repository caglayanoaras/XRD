''' This module is all about TiAlN'''

import pymatgen as mg
import subprocess
from pymatgen.io.lammps.outputs import parse_lammps_dumps,LammpsDump

from time import process_time

def create_lammps_structure_file(pymatgen_structure):
    ''' Attributes: pymatgen structure 
        Returns: a LAMMPS structure file
        '''
    with open('structure','w') as fdata:
        fdata.write('LAMMPS structure written for MC- Only TiAlN\n\n')

        fdata.write('{} atoms\n\n'.format(len(pymatgen_structure)))
        fdata.write('{} atom types\n\n'.format(3))

        fdata.write('{} {} xlo xhi\n'.format(0.0, pymatgen_structure.lattice.a))
        fdata.write('{} {} ylo yhi\n'.format(0.0, pymatgen_structure.lattice.b))
        fdata.write('{} {} zlo zhi\n'.format(0.0, pymatgen_structure.lattice.c))
        fdata.write('\n')

        fdata.write('Masses\n\n')
        fdata.write('{} {}\n'.format(1,26.981539))
        fdata.write('{} {}\n'.format(2,14.006700))
        fdata.write('{} {}\n\n'.format(3,47.867000))

        fdata.write('Atoms\n\n')

        for c,site in enumerate(pymatgen_structure):
            if site.species.formula == 'Ti1':
              fdata.write('{} {} {} {} {}\n'.format(c+1, 3,*site.coords))
            if site.species.formula == 'Al1':
              fdata.write('{} {} {} {} {}\n'.format(c+1, 1,*site.coords))      
            if site.species.formula == 'N1':
              fdata.write('{} {} {} {} {}\n'.format(c+1, 2,*site.coords))
              
def read_lammps_dump_file(filename = ''):
    ''' Attributes: path of a LAMMPS dump file
        Returns: pymatgen structure of final relaxed structure.
        '''
    dump = parse_lammps_dumps(filename)
    for i in dump:
        a = i
    lattice = a.box.to_lattice()
    df = a.data.sort_values(by=['id'])
    
    element_array = []
    for i in df['type'].tolist():
        if i ==1:
            element_array.append('Al')
        if i ==2:
            element_array.append('N')
        if i ==3:
            element_array.append('Ti')

    cord_np = df[['xs', 'ys', 'zs']].values      
    pymatgen_struct = mg.Structure(lattice, element_array, cord_np)
    return pymatgen_struct

    
def run_molecular_dynamics(pymatgen_structure):
    ''' Attributes: a pymatgen structure - perfect crystal preferably
        Returns: a pymatgen structure - relaxed version of input structure
        '''
    create_lammps_structure_file(pymatgen_structure)
    subprocess.run('lmp_serial.exe -in in.EVcurve', shell = True)
    relaxed_structure = read_lammps_dump_file('dump.atoms')

    return relaxed_structure

    
if __name__ == '__main__':
    
    file = 'SurfaceRandom.json'
    pymatgen_structure = mg.Structure.from_file(file)
    print(pymatgen_structure)
    print('---')
    t1_start = process_time()
    relaxed = run_molecular_dynamics(pymatgen_structure)
    t1_stop = process_time()
    relaxed.to(filename = 'relaxed_structure.json')
    print("Elapsed time:", t1_stop, t1_start)
    print("Elapsed time during the MD run in seconds:", 
                                         t1_stop-t1_start) 
    print(relaxed)

