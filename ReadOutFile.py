'''
    This module reads an out file and produces an excel file with Raw and
    Deconv sheets. 
'''
import XRDLib            as functions
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from pymatgen.core.spectrum import Spectrum

plt.style.use(['seaborn-paper', 'presentation'])

#-*-*-*-*-*-*-*-*-*-*-*-*-INPUT-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

fraction_inst        = 0.2     #Instrumental broadening fraction. Check alumina.py
wavelength           = 1.54056 # CuKa1 in Angstrom
deg_of_bck_poly      = 5      # Degree of background polynomial. Integer
name_of_exceloutput  = 'GADDS_960.xlsx'
name_of_outfileinput = 'outFiles//D5000//new_tube//TiAlN_560C_1mm_stripped.out'

#Enter expected peaks for smooth fitting
expected_peaks       = [37.52 , 43.58 , 63.29, 75.9 , 79.8, 95.6] #43.50 for 710

#Enter UVW parameters obtained from alumina.py module.
if 'D5000' in name_of_outfileinput:
    uvw = [ 0.00556691, -0.00089748,  0.00282308] #D5000
if 'GADDS' in name_of_outfileinput:
    uvw = [ 0.04586492, 0.0146818  ,  0.01856546] #GADDS
    
#-*-*-*-*-*-*-*-*-*-*-*-*-INPUT-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


#Read from out file
two_theta,intensity = functions.read_out_file(name_of_outfileinput)
index     = [c for c,i in enumerate(two_theta) if i>34 and i<100] # this bump is coming from sample holder.
intensity = [i for c,i in enumerate(intensity) if c in index]
two_theta = [i for c,i in enumerate(two_theta) if c in index]

#Normalize out file. Max intensity will be 100
spec = Spectrum(two_theta,intensity )
spec.normalize(mode = 'max',value = 100)

#This will fit on the scatter and return the lmfit result
result = functions.fit_experimental_data(spec.x,spec.y,expected_peaks,deg_of_bck_poly=deg_of_bck_poly)

# Show the fit
f, (a0, a1) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[3, 1]})
f.set_size_inches(8, 16)
a0.scatter(spec.x, spec.y, s = 8)
a0.plot(spec.x, result.best_fit, 'r-', linewidth=1.4)

residual = [a_i - b_i for a_i,b_i in zip(spec.y, result.best_fit)]
a1.plot(two_theta, residual, lw = 1.2)
a0.set_ylim([0,120])
a1.set_ylim([-10,10])
a0.set_xlim([35,100])
a1.set_xlim([35,100])

a0.grid( which='both', linestyle=':')
a0.set_xlabel('2\u03F4',fontsize = 14)
a0.set_ylabel('Intensity A.u.',fontsize = 14)

a1.set_ylabel('Residual',fontsize = 14)
a0.set_title('TiAlN powder - D5000 '.format(len(two_theta)),fontsize = 14)
a1.set_xticklabels([])

a0.xaxis.labelpad = 10
plt.show()
# Show the fit

import sys
sys.exit()
#keep the line profiles of sample as lists
fwhm,peak_two_thetas,amplitude,fraction, _ = functions.get_profile_data(two_theta, result)
d_hkl  = [wavelength/2/np.sin(np.pi/360*i) for i in peak_two_thetas]

al_fwhm     = [functions.caglioti(i,*uvw) for i in peak_two_thetas]
al_fraction = [fraction_inst] *len(al_fwhm)  # 0 for gaussian

#Deconvolution
deconv_fwhm            = []
deconv_fraction        = []
deconv_amplitude       = []
deconv_peak_two_thetas = []
def deconvolution():

    for i in range(6):
        _,deconv_popt = functions.get_deconvoluted_profile(fwhm,peak_two_thetas,amplitude,fraction,d_hkl,
                                                           al_fwhm,al_fraction, i, plot = False)
        deconv_fwhm.append(deconv_popt[1])
        deconv_fraction.append(deconv_popt[2])
        deconv_amplitude.append(deconv_popt[3])
        deconv_peak_two_thetas.append(deconv_popt[0])
        
def bruteDeconvolution():
    after_deconv   = []
    global name_of_exceloutput
    for c in range(len(peak_two_thetas)):
        x1x2 = functions.Deconvolver(fwhm[c], fraction[c], al_fwhm[c], al_fraction[c])
        after_deconv.append(x1x2)

    deconv_fwhm       = [row[0] for row in after_deconv]
    deconv_fraction   = [row[1] for row in after_deconv]
    deconv_amplitude  = [i for i in amplitude]
    deconv_peak_two_thetas = [i for i in peak_two_thetas]
    name_of_exceloutput = name_of_exceloutput[:-5] + '_bruteDeconv.xlsx'
    
deconvolution() # Select either deconvolution or bruteDeconvolution


# Print everyhing about measurement to an Excel file.
planes = [[(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)],
          [(2,0,0),(0,2,0),(0,0,2)],
          [(2,2,0),(2,0,2),(0,2,2),(2,-2,0),(2,0,-2),(0,2,-2)],
          [(3,1,1),(1,3,1),(1,1,3),(3,1,-1),(1,3,-1),(1,-1,3),(-3,1,-1),(1,-3,-1),(1,-1,-3),(3,-1,1),(-1,3,1),(-1,1,3)],
          [(2,2,2),(-2,2,2),(2,-2,2),(2,2,-2)],
          [(4,0,0),(0,4,0),(0,0,4)]]

raw_sample_sigma_deld_over_d    = [i*np.pi/180/4/np.tan(j*np.pi/360) for i,j in zip(fwhm,peak_two_thetas)]
deconv_sample_sigma_deld_over_d = [i*np.pi/180/4/np.tan(j*np.pi/360) for i,j in zip(deconv_fwhm,deconv_peak_two_thetas)]

Q        = [4*np.pi * np.sin(i*np.pi / 360) / wavelength for i in peak_two_thetas]
deconv_Q = [4*np.pi * np.sin(i*np.pi / 360) / wavelength for i in deconv_peak_two_thetas]

raw_sample_dict    = {'fwhm': fwhm, 'peak_two_thetas': peak_two_thetas, 'amplitude': amplitude, 'fraction': fraction,
                      'd_hkl':d_hkl, 'planes':planes, 'sigma_deltad_d':raw_sample_sigma_deld_over_d, 'Q': Q, 'sigma_of_d': [(i*j) for i,j in zip(raw_sample_sigma_deld_over_d, d_hkl)]}
deconv_sample_dict = {'fwhm': deconv_fwhm, 'peak_two_thetas': deconv_peak_two_thetas, 'amplitude': deconv_amplitude, 'fraction': deconv_fraction,
                      'd_hkl':d_hkl, 'planes':planes, 'sigma_deltad_d':deconv_sample_sigma_deld_over_d, 'Q': deconv_Q, 'sigma_of_d': [(i*j) for i,j in zip(deconv_sample_sigma_deld_over_d, d_hkl)]}

raw_sample_pd    = pd.DataFrame(raw_sample_dict)
deconv_sample_pd = pd.DataFrame(deconv_sample_dict)

writer = pd.ExcelWriter('ExcelSheets/' + name_of_exceloutput , engine='xlsxwriter')

raw_sample_pd.to_excel(writer, sheet_name='RawSample')
deconv_sample_pd.to_excel(writer, sheet_name='DeconvSample')

writer.save()

