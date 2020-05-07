import pymatgen as mg
import pandas   as pd
import numpy    as np
import matplotlib.pyplot as plt

from pymatgen.core.spectrum import Spectrum
from scipy.optimize         import curve_fit
from lmfit                  import Model
from lmfit.models           import PseudoVoigtModel, PolynomialModel
from numpy.fft              import rfft, rfftfreq, irfft

C_0 = 4* np.log(2)
C_1 = 4
pi = np.pi
indi = 101
wavelength = 1.54056

def Gaussian(x,theta_k,H_k):
    return np.sqrt(C_0)/(H_k * np.sqrt(pi)) * np.exp(-C_0*((x -  theta_k) / H_k)**2)

def Lorentzian(x,theta_k,H_k):
    return np.sqrt(C_1)/(H_k * pi) / (1 + C_1 * ((x - theta_k)/ H_k)**2)

def pseudo_voigt(x,theta_k,H_k,mu,amp):
    return (mu * Lorentzian(x,theta_k,H_k) + (1- mu) * Gaussian(x,theta_k,H_k))*amp

def read_out_file(path_of_file):
    ## reads out file and returns two_theta and intensity arrays
    two_theta = []
    intensity = []
    file = open(path_of_file, 'r')
    contents = file.readlines()
    for i in contents:
        if i[0] != '!':
            a = i.split(' ')
            two_theta.append(float(a[0]))
            intensity.append(float(a[1]))
    return two_theta,intensity

def read_xlsx_file(path_of_file):
    ## reads xlsx file and returns two_theta and intensity arrays
    df = pd.read_excel(path_of_file)
    two_theta = df.iloc[:, 0].astype(float)
    intensity = df.iloc[:, 1].astype(float)
    return two_theta,intensity

def fit_experimental_data(exp_x, exp_y,expected_peak_pos, deg_of_bck_poly = 5, maxfev = 25000):
    num_of_peaks = len(expected_peak_pos)
    mod = PolynomialModel(deg_of_bck_poly,prefix = 'poly_')    
    for c in range(num_of_peaks):
        mod = mod + PseudoVoigtModel(prefix = 'p{}_'.format(c))
    
    params = mod.make_params()

    center = 0
    sigma = 0
    amplitude = 0
    
    for param in params:
        if 'center' in param:
            params[param].set(value = expected_peak_pos[center])
            params[param].set(min = expected_peak_pos[center]-0.5)
            params[param].set(max = expected_peak_pos[center]+0.5)
            center += 1
        if 'poly' in param:
            if param == 'poly_c0':
                params[param].set(value = 0)
##                params[param].set(vary = False)
                
                params[param].set(min = -100)
                params[param].set(max = 100)
                continue
            if param == 'poly_c1':
                params[param].set(value = -1)
                params[param].set(min = -100)
                params[param].set(max = 100)
                continue
            params[param].set(value = 0)
##            params[param].set(min = 3e-1)
##            params[param].set(max = 3e-1)
        if 'sigma' in param:
            params[param].set(value = 0.5)
            params[param].set(min = 0.0001)
            params[param].set(max = 0.8)
            sigma += 1
        if 'amplitude' in param:
            params[param].set(value = 5.5)
            params[param].set(min = 0.0001)
            amplitude += 1
    result = mod.fit(np.asarray(exp_y),params, x=np.asarray(exp_x), fit_kws={'maxfev': maxfev})

    print(result.fit_report())
    return result

def get_profile_data(exp_x, result,excel = False):
# This func will take the result of experimental fit, clear the background, print an excel file with clean data and return clean intensity
    center = []
    amplitude = []
    fraction = []
    fwhm = []
    for name,value in zip(list(result.params.valuesdict().keys()), list(result.params.valuesdict().values())):
        if 'center' in name:
            center.append(value)
        if 'amplitude' in name:
            amplitude.append(value)
        if 'fraction' in name:
            fraction.append(value)
        if 'fwhm' in name:
            fwhm.append(value)
    clean_intensity = []
    for two_theta in exp_x:
        intensity = 0
        for c,cent in enumerate(center):
            intensity += pseudo_voigt(two_theta,cent,fwhm[c],fraction[c],amplitude[c])
        clean_intensity.append(intensity)
    spec = Spectrum(exp_x,clean_intensity)
    spec.normalize(mode = 'max', value = 100)
    if excel == True:
        dict_intense = {'two_theta' : spec.x , 'intensity': spec.y}
        dict_intense = pd.DataFrame(dict_intense)
        dict_intense.to_excel('data_400_exp_without_bck.xlsx')
    return fwhm, center,amplitude,fraction, spec.y

def caglioti(two_theta,u,v,w):
    theta = two_theta/2.0
    Hf = np.sqrt(u*(np.tan(2*np.pi / 360.0 * theta ))**2 + v * np.tan(2*np.pi / 360.0 * theta) + w)
    return Hf
def WH(two_theta,strain,L):
    theta = two_theta/2.0
    Hf = 180/np.pi*(4* strain * np.tan(np.pi / 180 * theta) + (1.54056/L / np.cos(np.pi / 180 * theta))) # fwhm in degrees
    return Hf
def linear(x,a,b):
    y = a*x + b
    return y
def Deconvolver(z1,z2,y1,y2):
# Inputs = 2 parameter of output signal -FWHM Fraction 
#           2 parameter of instrumental Broadening
# Outputs = 2 parameter of sample broadening
    popt = [0, z1, z2, 1]  # output with mean 0 and amplitude 1
    p02 =  [0, y1, y2, 1]  # instrumental boradening with mean 0 and amplitude 1
    
    xp = np.linspace(-4, 4, 4000)
    dx = xp[1] - xp[0]
    
    aa1 = np.linspace(0.1,1.5,200)
    aa2 = np.linspace(0.,1.,150)
    error = []
    a1a2= []
    for a2 in aa1:
        for a3 in aa2:
            a1a2.append([a2,a3])
            convolution = np.convolve(pseudo_voigt(xp ,0,a2,a3,1),pseudo_voigt(xp, *p02), mode="same")*dx
            poptt, _ = curve_fit(pseudo_voigt, xp, convolution,bounds=([0,0,0,1], [0.001, 1.5,1, 1.001]))
            error.append(sum([abs(a-b) for a,b in zip(popt,poptt)]))
            
    return a1a2[error.index(min(error))]
#### Added for figures
def fit_experimental_data_gauss(exp_x, exp_y,expected_peak_pos, deg_of_bck_poly = 5, maxfev = 25000):
    num_of_peaks = len(expected_peak_pos)
    mod = PolynomialModel(deg_of_bck_poly,prefix = 'poly_')    
    for c in range(num_of_peaks):
        mod = mod + PseudoVoigtModel(prefix = 'p{}_'.format(c))
    
    params = mod.make_params()

    center = 0
    sigma = 0
    amplitude = 0
    fraction = 0
    for param in params:
        if 'center' in param:
            params[param].set(value = expected_peak_pos[center])
            params[param].set(min = expected_peak_pos[center]-0.5)
            params[param].set(max = expected_peak_pos[center]+0.5)
            center += 1
        if 'poly' in param:
            if param == 'poly_c0':
                params[param].set(value = 50)
                params[param].set(min = -100)
                params[param].set(max = 100)
                continue
            if param == 'poly_c1':
                params[param].set(value = -1)
                params[param].set(min = -100)
                params[param].set(max = 100)
                continue
            params[param].set(value = 0)
##            params[param].set(min = 3e-1)
##            params[param].set(max = 3e-1)
        if 'sigma' in param:
            params[param].set(value = 0.5)
            params[param].set(min = 0.0001)
            params[param].set(max = 0.8)
            sigma += 1
        if 'amplitude' in param:
            params[param].set(value = 5.5)
            params[param].set(min = 0.0001)
            amplitude += 1
        if 'fraction' in param:
            params[param].set(value = 0.0)
            params[param].set(min = 0.000)
            params[param].set(max = 0.000001)
            fraction += 1
    result = mod.fit(np.asarray(exp_y),params, x=np.asarray(exp_x), fit_kws={'maxfev': maxfev})

    print(result.fit_report())
    return result
def fit_experimental_data_lorentz(exp_x, exp_y,expected_peak_pos, deg_of_bck_poly = 5, maxfev = 25000):
    num_of_peaks = len(expected_peak_pos)
    mod = PolynomialModel(deg_of_bck_poly,prefix = 'poly_')    
    for c in range(num_of_peaks):
        mod = mod + PseudoVoigtModel(prefix = 'p{}_'.format(c))
    
    params = mod.make_params()

    center = 0
    sigma = 0
    amplitude = 0
    fraction = 0
    for param in params:
        if 'center' in param:
            params[param].set(value = expected_peak_pos[center])
            params[param].set(min = expected_peak_pos[center]-0.5)
            params[param].set(max = expected_peak_pos[center]+0.5)
            center += 1
        if 'poly' in param:
            if param == 'poly_c0':
                params[param].set(value = 50)
                params[param].set(min = -100)
                params[param].set(max = 100)
                continue
            if param == 'poly_c1':
                params[param].set(value = -1)
                params[param].set(min = -100)
                params[param].set(max = 100)
                continue
            params[param].set(value = 0)
##            params[param].set(min = 3e-1)
##            params[param].set(max = 3e-1)
        if 'sigma' in param:
            params[param].set(value = 0.5)
            params[param].set(min = 0.0001)
            params[param].set(max = 0.8)
            sigma += 1
        if 'amplitude' in param:
            params[param].set(value = 5.5)
            params[param].set(min = 0.0001)
            amplitude += 1
        if 'fraction' in param:
            params[param].set(value = 1.0)
            params[param].set(min = 1.000)
            params[param].set(max = 1.000001)
            fraction += 1
    result = mod.fit(np.asarray(exp_y),params, x=np.asarray(exp_x), fit_kws={'maxfev': maxfev})

    print(result.fit_report())
    return result

### This part is added for deconvolution of instrumental broadening
def profile_func(measurement,LP, pointnum):

    global a
    a = 2*np.pi
    h3 = 2*a/wavelength*np.sin(measurement['peak_two_thetas'][LP]/360*np.pi)
    theta1 = np.arcsin((h3-0.5)*wavelength/2/a)*360/np.pi
    theta2 = np.arcsin((h3+0.5)*wavelength/2/a)*360/np.pi
    
    x = np.linspace(theta1,theta2, num = pointnum)
    y = pseudo_voigt(x,measurement['peak_two_thetas'][LP],
                      measurement['fwhm'][LP],
                      measurement['fraction'][LP],
                      measurement['amplitude'][LP])
        
    profile2theta = Spectrum(x,y)
    return profile2theta

def profile_func_2theta_to_h3(spectrum2theta):
    '''
    Args: spectrum2theta > pymatgen spectrum of 2 theta
    Returns: pymatgen spectrum of h3
    '''
    cos = np.cos(spectrum2theta.x *np.pi/360)
    y = np.divide(spectrum2theta.y, cos) * wavelength/a
    x = 2*a/wavelength* np.sin(spectrum2theta.x *np.pi/360)
    profileh = Spectrum(x,y)
    popt,_ = curve_fit(pseudo_voigt, profileh.x, profileh.y,bounds = ([profileh.x.min(), 0, 0, 0,],
                                                                      [profileh.x.max(), 0.5, 1, 240 ]))
    return profileh,popt
def profile_func_h3_to_2theta(profileh):

    x = 360/np.pi*np.arcsin(profileh.x*wavelength/2/a )
    
    cos = np.cos(x *np.pi/360)
    y = np.multiply(profileh.y, cos) / wavelength*a
    
    #
    spectrum2theta = Spectrum(x,y)
    popt,_ = curve_fit(pseudo_voigt, spectrum2theta.x, spectrum2theta.y,bounds = ([spectrum2theta.x.min(), 0.2, 0, 0,],
                                                                      [spectrum2theta.x.max(), 2, 1, 240 ]))
    return spectrum2theta,popt

def get_deconvoluted_profile(fwhm,peak_two_thetas,amplitude,fraction,d_hkl,al_fwhm,al_fraction, index, plot = True):
    meas_dict = {'fwhm':fwhm, 'peak_two_thetas':peak_two_thetas, 'amplitude':amplitude, 'fraction':fraction, 'd_hkl':d_hkl }
    al_dict   = {'fwhm':al_fwhm, 'peak_two_thetas':peak_two_thetas, 'amplitude':amplitude, 'fraction':al_fraction, 'd_hkl':d_hkl }
    measurement = pd.DataFrame(meas_dict)
    instrument  = pd.DataFrame(al_dict)
    mesh       = 40001
    N          = 40001
    
    profile2theta            = profile_func(measurement,index, mesh)
    profile2theta_instrument = profile_func(instrument,index, mesh)
    
    profileh,popt                       = profile_func_2theta_to_h3(profile2theta)
    profileh_instrument,popt_instrument = profile_func_2theta_to_h3(profile2theta_instrument)

    
    x= np.linspace(profileh.x.min(),profileh.x.max(),num=N)
    h = pseudo_voigt(x,*popt)
    h_inst = pseudo_voigt(x,*popt_instrument)


    H = rfft(h)
    H_INST = rfft(h_inst)
    frequency = rfftfreq(N,round(x[1]-x[0],6))
    
    ind = np.arange(0,indi)
    freq = frequency[ind]
    normalizer = 1#1/abs(H[0])
    coefs = np.abs(H[ind])*normalizer
    coefs_inst = np.abs(H_INST[ind])*normalizer
    coefs_deconv = coefs/coefs_inst *coefs[0]
    deconv_q = np.roll(irfft(coefs_deconv, n =N),N//2)
    deconvq_spectrum = Spectrum(x,deconv_q)
    if plot == True:
##        plt.scatter(freq, np.imag(H[ind]*normalizer), label = 'imaginary part')
        plt.scatter(freq, coefs, label = 'measured profile' )
        plt.scatter(freq, coefs_inst, label = 'instrumental profile' )
        plt.scatter(freq, coefs_deconv, label = 'deconvoluted profile' )
        plt.title('plane index: {}'.format(index))
        plt.xlabel('n')
        plt.ylabel('An')
        plt.legend()
        plt.grid()
        plt.show()
        
##        f, (ax1,ax2) = plt.subplots(1,2)
##        ax1.scatter(profile2theta.x, profile2theta.y, s = 6)
##        ax1.set_xlabel('2\u03B8')
##        ax2.scatter(profileh.x, profileh.y, s = 6)
##        ax2.plot(x, pseudo_voigt(x, *popt), label = 'measured profile',alpha = 0.8)
##        ax2.plot(x, deconv_q, label='deconvoluted profile',alpha = 0.8)
##        ax2.set_xlabel('h')
##        ax1.grid()
##        ax2.grid()
##        ax2.legend()
##        plt.show()
    
    return profile_func_h3_to_2theta(deconvq_spectrum)
