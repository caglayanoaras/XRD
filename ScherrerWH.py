'''
    This module reads an excel file and produces Shrerrer and WH results.
'''
import XRDLib            as functions
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mlt
from pymatgen.core.spectrum import Spectrum
from sklearn.preprocessing  import PolynomialFeatures
from sklearn.linear_model   import LinearRegression
from scipy.optimize         import curve_fit
from os                     import listdir
from os.path                import isfile, join

plt.style.use(['seaborn-paper', 'presentation'])
plt.rc('font', family='serif')

wavelength         = 1.54056

excel_path         = 'excelSheets/GADDS_1060.xlsx'
measurement_deconv = pd.read_excel(excel_path, sheet_name='DeconvSample', index_col=0)
measurement        = pd.read_excel(excel_path, sheet_name='RawSample', index_col=0)

# Williamson-Hall analysis for deconvoluted sample. Figure also includes raw sample fwhms
xforWH_measured = [4*np.sin(i*np.pi/360) for i in measurement['peak_two_thetas']]
xforWH_material = [4*np.sin(i*np.pi/360) for i in measurement['peak_two_thetas']] # after deconvolution

yforWH_measured = [j*np.pi/180*np.cos(i*np.pi/360) for i,j in zip(measurement['peak_two_thetas'],measurement['fwhm'])]
yforWH_material = [j*np.pi/180*np.cos(i*np.pi/360) for i,j in zip(measurement_deconv['peak_two_thetas'],measurement_deconv['fwhm'])] # after deconvolution

# Show the fit
plt.figure(figsize=(12,8))
#plt.title('Williamson-Hall Fit \n')
#plt.scatter(xforWH_measured, yforWH_measured, color = 'k') #,label = 'Measured FWHM'
plt.scatter(xforWH_material, yforWH_material, color = 'r' ) #,label = 'Material FWHM'
plt.xlabel('4sin\u03F4')
plt.ylabel('FWHM(rad)cos\u03F4')
plt.grid(linestyle='--')
plt.xlim(0,4)
plt.ylim(0,None)
popt, pcov = curve_fit(functions.linear, xforWH_material, yforWH_material,bounds=([0,0], [3., 20000000.]))
xmesh = np.linspace(0,4,100)
plt.plot(xmesh, functions.linear(xmesh, *popt), 'k-',label = 'strain + crys. size fit')#,label = 'WH fit'

# Print Results
print('Williamson Hall results:')
print('Strain (\u03B5)          : {} %'.format(popt[0]*100))
print('Crystallite Size (L): {} Angstrom - K is taken 0.94'.format(0.94*wavelength/popt[1]))
print('-*-*-*-*-*-*-')
popt, pcov = curve_fit(functions.linear, xforWH_material, yforWH_material,bounds=([0,0], [3., 1e-16]))
plt.plot(xmesh, functions.linear(xmesh, *popt), 'b--',label = 'Only strain fit')
print('Strain (\u03B5) - only strain          : {} %'.format(popt[0]*100))
print('Crystallite Size (L) - only strain : {} Angstrom - K is taken 0.94'.format(0.94*wavelength/popt[1]))
print('-*-*-*-*-*-*-')
popt, pcov = curve_fit(functions.linear, xforWH_material, yforWH_material,bounds=([0,0], [1e-16, 20000000.]))
plt.plot(xmesh, functions.linear(xmesh, *popt), 'g--',label = 'Only crys. size fit')
print('Strain (\u03B5) - only crystallite          : {} %'.format(popt[0]*100))
print('Crystallite Size (L) - only crystallite : {} Angstrom - K is taken 0.94'.format(0.94*wavelength/popt[1]))
print('-*-*-*-*-*-*-')
# Print Results

plt.legend()
plt.show()
# Show the fit


# Scherrer equation only for 111 and 200 planes
def scherrer(fwhm, theta):
    K = 0.94
    L = K*wavelength/fwhm/np.cos(theta)
    return L

L111_deconv = scherrer(measurement_deconv['fwhm'][0]*np.pi/180, measurement_deconv['peak_two_thetas'][0]*np.pi/360)
L200_deconv = scherrer(measurement_deconv['fwhm'][1]*np.pi/180, measurement_deconv['peak_two_thetas'][1]*np.pi/360)
L111        = scherrer(measurement['fwhm'][0]*np.pi/180, measurement['peak_two_thetas'][0]*np.pi/360)
L200        = scherrer(measurement['fwhm'][1]*np.pi/180, measurement['peak_two_thetas'][1]*np.pi/360)

print('Scherrer equation on 111 profile -no deconv.  -: {}'.format(L111))
print('Scherrer equation on 111 profile -with deconv.-: {}'.format(L111_deconv))

print('Scherrer equation on 200 profile -no deconv.  -: {}'.format(L200))
print('Scherrer equation on 200 profile -with deconv.-: {}'.format(L200_deconv))




