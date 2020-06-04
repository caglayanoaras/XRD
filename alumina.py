# This script will read alumina out files and get caglioti equation parameters for the profile fit
# Does not return any file.

import XRDLib            as functions
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mlt

from pymatgen.core.spectrum import Spectrum
from sklearn.preprocessing  import PolynomialFeatures
from sklearn.linear_model   import LinearRegression
from scipy.optimize         import curve_fit

wavelength = 1.54056 # CuKa1 in Angstrom

plt.style.use(['seaborn-paper', 'presentation'])

# Reading out files
two_theta_D5000,intensity_D5000 = functions.read_out_file('outFiles//D5000//1mm//Al2O3_2h30min_stripped.out')
two_theta_GADDS,intensity_GADDS = functions.read_out_file('outFiles//GADDS//alumina_15min.out')


# In order make a better fit some intensity values are changed to 0. After that intensities are normalized that maximum value become 100.
indexes_to_make_zero = []
for c,i in enumerate(two_theta_D5000):       
    if 41<i<42 or 61<i<62 or 73<i<88.5 or 0<i<25 or 89.5<i<94 or (intensity_D5000[c]<0.45):
        indexes_to_make_zero.append(c)    
for c in indexes_to_make_zero:
    intensity_D5000[c] = 0
spec_D5000 = Spectrum(two_theta_D5000,intensity_D5000 )
spec_D5000.normalize(mode = 'max',value = 100)


#Expected peaks are selected for the fit.
expected_peaks_D5000 = [25.6,35.1,37.79,43.4,52.6,57.6,66.56,68.25,89.1,95.34, 101.14]#,61.35


#Fit is done on scatter data. 
result_D5000 = functions.fit_experimental_data(spec_D5000.x,spec_D5000.y,expected_peaks_D5000,deg_of_bck_poly = 0)


# Show the fit
f1, (a0, a1) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[3, 1]})
f1.set_size_inches(8, 16)
a0.scatter(spec_D5000.x, spec_D5000.y, s = 8)
a0.plot(spec_D5000.x, result_D5000.best_fit, 'r-', linewidth=1.4)

residual_D5000 = [a_i - b_i for a_i,b_i in zip(spec_D5000.y, result_D5000.best_fit)]
a1.plot(spec_D5000.x, residual_D5000, lw = 1.4)
a0.set_ylim([0,120])
a1.set_ylim([-10,10])
a0.set_xlim([20,105])
a1.set_xlim([20,105])

a0.grid( which='both', linestyle=':')
a0.set_xlabel('2\u03F4',fontsize = 14)
a0.set_ylabel('Normalized Intensity',fontsize = 14)

a1.set_ylabel('Residual',fontsize = 14)
a0.set_title('Fit of Alumina - D5000  ',fontsize = 14)
a1.set_xticklabels([])

a0.xaxis.labelpad = 10
plt.show()
# Show the fit ended.


#Keep the line profiles of sample as lists
fwhm_D5000,peak_two_thetas_D5000,amplitude_D5000,fraction_D5000, _ = functions.get_profile_data(two_theta_D5000, result_D5000)


#Show FWHM versus 2theta of alumina '''
uvw_D5000, _ = curve_fit(functions.caglioti, peak_two_thetas_D5000, fwhm_D5000, maxfev=2000)
f2, ax1 = plt.subplots()
ax1.scatter(peak_two_thetas_D5000 , fwhm_D5000 )
ax1.plot(spec_D5000.x, [functions.caglioti(i, *uvw_D5000) for i in spec_D5000.x], 'r-')
ax1.set_title('D5000 Measurement - Alumina')
ax1.set_ylabel('FWHM (\u00B0)')
ax1.set_xlabel('2\u03F4 (\u00B0)')
ax1.grid( which='both', linestyle=':')
plt.show()
print('Considering 10 peaks of alumina, UVW parameters for D5000:')
print(uvw_D5000)


## D5000 is done. Same procedure is taken for GADDS.

# In order make a better fit some intensity values are changed to 0. After that intensities are normalized that maximum value become 100.
indexes_to_make_zero = []
for c,i in enumerate(two_theta_GADDS):       
    if 41<i<42 or 61<i<62 or 73<i<88.5 or 89.9<i<94:#or 0<i<24.5 
        indexes_to_make_zero.append(c)    
for c in indexes_to_make_zero:
    intensity_GADDS[c] = 2 #background is around 7
spec_GADDS = Spectrum(two_theta_GADDS,intensity_GADDS )
spec_GADDS.normalize(mode = 'max',value = 100)


#Expected peaks are selected for the fit.
expected_peaks_GADDS = [25.6,35.1,37.75,43.4,52.6,57.6,66.56,68.25,89.07,95.34, 101.14]#,61.35


#Fit is done on scatter data. 
result_GADDS = functions.fit_experimental_data(spec_GADDS.x,spec_GADDS.y,expected_peaks_GADDS,deg_of_bck_poly = 2)


# Show the fit
f3, (a0, a1) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[3, 1]})
f3.set_size_inches(8, 16)
a0.scatter(spec_GADDS.x, spec_GADDS.y, s = 8)
a0.plot(spec_GADDS.x, result_GADDS.best_fit, 'r-', linewidth=1.4)

residual_GADDS = [a_i - b_i for a_i,b_i in zip(spec_GADDS.y, result_GADDS.best_fit)]
a1.plot(spec_GADDS.x, residual_GADDS, lw = 1.4)
a0.set_ylim([0,120])
a1.set_ylim([-10,10])
a0.set_xlim([20,105])
a1.set_xlim([20,105])

a0.grid( which='both', linestyle=':')
a0.set_xlabel('2\u03F4',fontsize = 14)
a0.set_ylabel('Normalized Intensity',fontsize = 14)

a1.set_ylabel('Residual',fontsize = 14)
a0.set_title('Fit of Alumina - GADDS  ',fontsize = 14)
a1.set_xticklabels([])

a0.xaxis.labelpad = 10
plt.show()
# Show the fit ended.


#Keep the line profiles of sample as lists
fwhm_GADDS,peak_two_thetas_GADDS,amplitude_GADDS,fraction_GADDS, _ = functions.get_profile_data(two_theta_GADDS, result_GADDS)


#Show FWHM versus 2theta of alumina '''
uvw_GADDS, _ = curve_fit(functions.caglioti, peak_two_thetas_GADDS, fwhm_GADDS, maxfev=2000)
f4, ax1 = plt.subplots()
ax1.scatter(peak_two_thetas_GADDS , fwhm_GADDS )
ax1.plot(spec_GADDS.x, [functions.caglioti(i, *uvw_GADDS) for i in spec_GADDS.x], 'r-')
ax1.set_title('GADDS Measurement - Alumina')
ax1.set_ylabel('FWHM (\u00B0)')
ax1.set_xlabel('2\u03F4 (\u00B0)')
ax1.grid( which='both', linestyle=':')
plt.show()
print('Considering 10 peaks of alumina, UVW parameters for GADDS:')
print(uvw_GADDS)


## D5000 and GADDS are done know you have uvw parameters. Final plot comparing both.
#Show FWHM versus 2theta of alumina '''
f5, ax1 = plt.subplots()

ax1.scatter(peak_two_thetas_GADDS , fwhm_GADDS , color = 'b', label = 'GADDS')
ax1.plot(spec_GADDS.x, [functions.caglioti(i, *uvw_GADDS) for i in spec_GADDS.x], 'r-', label = 'Caglioti Fit')

ax1.scatter(peak_two_thetas_D5000 , fwhm_D5000, color = 'k', label = 'D5000' )
ax1.plot(spec_D5000.x, [functions.caglioti(i, *uvw_D5000) for i in spec_D5000.x], 'r-')

ax1.set_title('GADDS vs D5000 Measurement - Alumina')
ax1.set_ylabel('FWHM (\u00B0)')
ax1.set_xlabel('2\u03F4 (\u00B0)')
ax1.grid( which='both', linestyle=':')
plt.legend()
plt.show()

