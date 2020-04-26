# This piece of code is slightly different from pymatgen.

import numpy as np
import collections

WAVELENGTHS = {
    "CuKa": 1.54184,
    "CuKa2": 1.54439,
    "CuKa1": 1.54056,
}
import json

with open('atomic_scattering_params_ict.json') as f:
    ATOMIC_SCATTERING_PARAMETERS_ict = json.load(f)
with open('atomic_scattering_params_dft.json') as f:
    ATOMIC_SCATTERING_PARAMETERS_dft = json.load(f)
    
from pymatgen.core.spectrum import Spectrum

class DiffractionPattern(Spectrum):
    XLABEL = "$2\\Theta$"
    YLABEL =  "Intensity"

    def __init__(self,x,y,hkls,d_hkls):
        super().__init__(x,y,hkls,d_hkls)
        self.hkls = hkls
        self.d_hkls = d_hkls

from abc import ABC, abstractmethod

class AbstractDiffractionPatternCalculator(ABC):
    
    TWO_THETA_TOL  = 0.00001
    SCALED_INTENSITY_TOL = 0.001
    
    @abstractmethod
    def get_pattern(self,structure,scaled = True, two_theta_range = (0,90)):
        pass
    
    def get_unique_families(hkls):
        def is_perm(hkl1, hkl2):
            h1 = np.abs(hkl1)
            h2 = np.abs(hkl2)
            return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

        unique = collections.defaultdict(list)
        for hkl1 in hkls:
            found = False
            for hkl2 in unique.keys():
                if is_perm(hkl1, hkl2):
                    found = True
                    unique[hkl2].append(hkl1)
                    break
            if not found:
                unique[hkl1].append(hkl1)

        pretty_unique = {}
        for k, v in unique.items():
            pretty_unique[sorted(v)[-1]] = len(v)

        return pretty_unique
    
from math import sin, cos, asin, pi, degrees, radians

class XRDCalculator(AbstractDiffractionPatternCalculator):

    
    def __init__(self,wavelength = 'CuKa', debye_waller_factors = None, typ = 'ict'):
        self.radiation = wavelength
        self.wavelength = WAVELENGTHS[wavelength]
        self.debye_waller_factors = debye_waller_factors or {}
        self.SCALED_INTENSITY_TOL = 0.001
        self.AVAILABLE_RADIATION =   tuple(WAVELENGTHS.keys())
        self.typ = typ
    def get_pattern(self,structure,scaled = True, two_theta_range = (0,90)):
        wavelength = self.wavelength
        latt = structure.lattice
        is_hex = latt.is_hexagonal()
        if two_theta_range is None:
            min_r, max_r = (0, 2/ wavelength)
        else:
            min_r, max_r = [2 * sin(radians(t / 2)) / wavelength for t in two_theta_range]
        recip_lattice = latt.reciprocal_lattice_crystallographic
        recip_pts =  recip_lattice.get_points_in_sphere([[0, 0, 0]], [0, 0, 0], max_r)
        if min_r:
            recip_pts = [pt for pt in recip_pts if pt[1] >= min_r]
##        print(max_r)
##        for i in range(0,len(recip_pts)):
##            print(recip_pts[i])
        
        zs = []
        coeffs = []
        fcoords = []
        occus = []
        dwfactors = []
        for site in structure:
            for sp, occu in site.species_and_occu.items():
                zs.append(sp.Z)
                try:
                    if self.typ == 'ict':
                        c = ATOMIC_SCATTERING_PARAMETERS_ict[sp.symbol]
                    if self.typ == 'dft':
                        c = ATOMIC_SCATTERING_PARAMETERS_dft[sp.symbol]
                except KeyError:
                    raise ValueError('Unable to calculate XRD pattern')
                coeffs.append(c)
                dwfactors.append(self.debye_waller_factors.get(sp.symbol, 0))
                fcoords.append(site.frac_coords)
                occus.append(occu)
        zs = np.array(zs)
        coeffs = np.array(coeffs)
        fcoords = np.array(fcoords)
        occus = np.array(occus)
        dwfactors = np.array(dwfactors)
        peaks = {}
        two_thetas = []

        for hkl,g_hkl, ind, _ in sorted(recip_pts, key=lambda i: (i[1], -i[0][0], -i[0][1], -i[0][2])):
            hkl = [int(round(i)) for i in hkl]
            if g_hkl != 0:
                d_hkl = 1/g_hkl
                theta = asin(wavelength * g_hkl /2)
                s = g_hkl /2
                s2 = s**2

                g_dot_r = np.dot(fcoords, np.transpose([hkl])).T[0]
                fs = np.sum(coeffs[:,:,0]*np.exp(-coeffs[:,:,1]*s2),axis=1)

                dw_correction = np.exp(-dwfactors * s2)

                f_hkl = np.sum(fs*occus*np.exp(2j * pi *g_dot_r)*dw_correction)

                lorentz_factor = (1 + cos(2*theta)**2) /  (sin(theta)**2 * cos(theta))

                i_hkl = (f_hkl * f_hkl.conjugate()).real

                two_theta = degrees(2*theta)

                if is_hex:
                    hkl = (hkl[0], hkl[1], - hkl[0] - hkl[1], hkl[2])

                ind = np.where(np.abs(np.subtract(two_thetas, two_theta)) < 0.00001)
                ''' two_theta_tol'''
                
                if len(ind[0]) > 0:
                    peaks[two_thetas[ind[0][0]]][0] += i_hkl * lorentz_factor
                    peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
                else:
                    peaks[two_theta] = [i_hkl * lorentz_factor, [tuple(hkl)],d_hkl]
                    two_thetas.append(two_theta)
                    
        max_intensity = max([v[0] for v in peaks.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        for k in sorted(peaks.keys()):
            v = peaks[k]
            fam = AbstractDiffractionPatternCalculator.get_unique_families(v[1])
            if v[0] / max_intensity * 100 > self.SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v[0])
                hkls.append([{"hkl": hkl, "multiplicity": mult} for hkl, mult in fam.items()])

                d_hkls.append(v[2])

        xrd = DiffractionPattern(x, y, hkls, d_hkls)
        if scaled:
            xrd.normalize(mode = 'max', value = 100)
        return xrd
                
            
