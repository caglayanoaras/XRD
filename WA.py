'''
    This module reads an excel file and does warren-averbach analysis.
'''

import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from XRDLib                 import pseudo_voigt
from pymatgen.core.spectrum import Spectrum
from numpy.fft              import rfft, rfftfreq
from scipy.optimize         import curve_fit
from sklearn.linear_model   import LinearRegression
from scipy.optimize         import least_squares

pi = np.pi
e  = np.e 
plt.style.use(['seaborn-paper', 'presentation'])
plt.rc('font', family='serif')

DECONV      = pd.read_excel('excelSheets/GADDS_asdep.xlsx',sheet_name='DeconvSample')
NOT_DECONV  = pd.read_excel('excelSheets/GADDS_asdep.xlsx',sheet_name='RawSample')

measurement = DECONV

wavelength = 1.54056
indi = 401
 
def profile_func(measurement,LP, pointnum):
    '''
    Args: measurement > pandas df of measurement
          pointnum    > number of datapoints in x axis
          LP          > index of desired profile
    Returns: pymatgen spectrum where x is 2theta y is intensity 
    '''
    global a
    if LP in [0,1,2,3]:
        a = measurement['d_hkl'][LP]
    if LP in [4]:
        a = measurement['d_hkl'][0]
    if LP in [5]:
        a = measurement['d_hkl'][1]
    h3 = 2*a/wavelength*np.sin(measurement['peak_two_thetas'][LP]/360*np.pi)
    theta1 = np.arcsin((h3-0.5)*wavelength/2/a)*360/np.pi
    theta2 = np.arcsin((h3+0.5)*wavelength/2/a)*360/np.pi
    
    x = np.linspace(theta1,theta2, num = pointnum)
    y = pseudo_voigt(x,measurement['peak_two_thetas'][LP],
                      measurement['fwhm'][LP],
                      measurement['fraction'][LP],
                      100)
##                      measurement['amplitude'][LP])
                     
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

def get_fourier_coefs(measurement, index, plot = True):
    '''
    Args: measurement > pymatgen spectrum of 2 theta
          index       > index of profile 
    Returns: fourier coefficients A_n of profile
    '''
    mesh       = 5001
    N          = 5001
    
    profile2theta = profile_func(measurement,index, mesh)
    profileh,popt = profile_func_2theta_to_h3(profile2theta)
    print('index {}'.format(index))
    print('popt {}'.format(popt))
    x= np.linspace(profileh.x.min(),profileh.x.max(),num=N)
    h = pseudo_voigt(x,*popt)

    H = rfft(h)
    frequency = rfftfreq(N,round(x[1]-x[0],6))
    ind = np.arange(0,indi)
    freq = frequency[ind]
    normalizer = 1/abs(H[0])
    coefs = np.abs(H[ind])*normalizer
    if plot == True:
        plt.scatter(freq, np.imag(H[ind]*normalizer))
        plt.scatter(freq, coefs )
        plt.xlabel('n')
        plt.ylabel('An')
        plt.grid()
        plt.show()
        
        f, (ax1,ax2) = plt.subplots(1,2)
        ax1.scatter(profile2theta.x, profile2theta.y, s = 6)
        ax1.set_xlabel('2\u03B8')
        ax1.set_ylabel('A.u.')
        ax2.scatter(profileh.x, profileh.y, s = 6)
        ax2.plot(profileh.x, pseudo_voigt(profileh.x, *popt))
        ax2.set_xlabel('h3')
        ax1.grid()
        ax2.grid()
        plt.show()
    return freq[:indi],coefs[:indi]

def solve_log_normal(larea, lvolume):
    Larea = larea
    Lvolume = lvolume
    def f(variables):
        D_0, sigma = variables

        if sigma > 0:
            first_eq = 2/3*D_0*np.exp(5/2*(np.log(sigma)*np.log(sigma))) - Larea
            second_eq = 3/4*D_0*np.exp(7/2*(np.log(sigma)*np.log(sigma))) - Lvolume
            return [first_eq, second_eq]
        return 0
    try:
        solution = least_squares(f, (15, 1), bounds = ((0, 0), (50, 5)),ftol = 1e-10,xtol = 1e-10,gtol = 1e-10)
##        print(solution)
    except:
        print('Cannot fit')
        solution = [0,0]
    return solution

def log_normal(D_0,sigma,x):
    return 1/(np.sqrt(2*np.pi)*x*np.log(sigma)) * np.exp(-0.5 * (np.log(x/D_0)/np.log(sigma))**2)

def WA(plane_dir ):
    
    if plane_dir == 111:
        freq1, coefs1 = get_fourier_coefs(measurement, 0, plot = True)
        freq2, coefs2 = get_fourier_coefs(measurement, 4, plot = False)
    if plane_dir == 200:
        freq1, coefs1 = get_fourier_coefs(measurement, 1, plot = True)    
        freq2, coefs2 = get_fourier_coefs(measurement, 5, plot = False)
    
    coefslog1 = np.log(coefs1)
    coefslog2 = np.log(coefs2)
## This plot will plot the fourier coefficients of a profile
    plt.figure(figsize=(10,8))
##    plt.title('Fourier coefficients for {'+ '{}'.format(plane_dir)+ '}')
    msZn     = []
    msStrain = [0]
    count = 0
    for i,j in zip(coefslog1, coefslog2):
        if count <6:
            plt.plot([1, 4], [i, j], 'ko--')
            count+=1
        msZn.append((i-j)/3/2/(np.pi**2))
    for n,Zn in enumerate(msZn[1:]):
        msStrain.append(Zn/(n+1)/(n+1))
        
    msZn[0]   = 0.0 # gets negative value yet so close to zero
    rmsZn     = [np.sqrt(i) for i in msZn]
    rmsStrain = [np.sqrt(i) for i in msStrain]
    plt.ylabel('$\ln(A_{n})$')
    plt.xlabel('$l^2$')
    
    if plane_dir == 111:
        pos = [0.002, -0.018, -0.05, -0.08, -0.11, -0.14]
    if plane_dir == 200:
        pos = [0.002, -0.028, -0.08, -0.125, -0.17, -0.22]
        
    for i in range(6):
        plt.text(3.6, pos[i],'n={}'.format(i), fontsize=18)
    plt.grid()
    plt.show()
    print('First 10 avg. strain in' + '{}'.format(plane_dir) + 'direction:\n {}'.format(rmsStrain[:10])  )
    print('First 10 avg. Zn in' + '{}'.format(plane_dir) + 'direction:\n {}'.format(rmsZn[:10])  )
    print('*-*-*-*')

    Asize = np.exp(np.array(msZn)*2*1*(np.pi**2) +np.log(coefs1))

## This plot will plot the size coefficients of a profile    
    plt.figure(figsize=(8,8))
##    plt.title('Size coefficients for {'+ '{}'.format(plane_dir)+ '}')
    plt.xlabel('$n$')
    plt.ylabel('$A_{s}$')
    plt.scatter(range(indi),Asize,color = 'k')   
    reg = LinearRegression().fit( np.arange(5).reshape(-1,1),Asize[:5])   
    dummyx = np.linspace(0,-reg.intercept_/reg.coef_,num=50)
    dummyy = reg.coef_[0]*dummyx + reg.intercept_
    plt.plot(dummyx, dummyy, 'r--')
    plt.scatter(-reg.intercept_/reg.coef_,0,s=35)
    plt.grid()
    plt.show()
    
    return rmsStrain

if __name__ == '__main__':
       
##    strain111 = WA(111)
    strain200 = WA(200)
    
## This part is about strain versus L plot
    plt.figure(figsize = (14,6))
##    plt.title('Strain over length L')
    plt.xlabel('L (Å)')
    plt.ylabel('Strain (%)')
    plt.grid()
##    plt.axhline(y=0.39, color ='r', linestyle= '--',label='WH')
    LL = [i*a for i in range(15)]
    plt.xticks(LL)
##    plt.plot(LL[1:],[i*100 for i in strain111[1:15]],'o-',label='111')
    plt.plot(LL[1:],[i*100 for i in strain200[1:15]],'o-')#,label='200')

    plt.legend()
    plt.show()



    


##    Larea    = round(float(-reg.intercept_/reg.coef_), 4)*a/10
##    Lvolume  = round(np.trapz(Asize)*2, 4)*a/10
##
##    print('<L>area   ' + '{}'.format(plane_dir) + ' direction: {}'.format(Larea))
##    print('<L>volume ' + '{}'.format(plane_dir) + ' direction: {}'.format(Lvolume))      
##    print('*-*-*-*')

    
## This part is about crystallite size distribution
####    D_0, sigma = solve_log_normal(Larea, Lvolume).x
####    print('D_0(median)    , ' + '{}'.format(plane_dir) + ' plane : {}'.format(round(D_0,4)))
####    print('Sigma(variance), ' + '{}'.format(plane_dir) + ' plane : {}'.format(round(sigma,4)))
####
####    print('*-*-*-*')
####    
####    print('<D>area   ' + '{}'.format(plane_dir) + ' direction: {}'.format(D_0*np.exp(5/2*np.log(sigma)*np.log(sigma))))
####    print('<D>volume ' + '{}'.format(plane_dir) + ' direction: {}'.format(D_0*np.exp(7/2*np.log(sigma)*np.log(sigma))))
####    print('*-*-*-*')
####
####    xforlognormal    = np.linspace(0.0001, 40, num = 1500)
####    yforlognormal = log_normal(D_0, sigma, xforlognormal)
####    
####    plt.figure(figsize=(14,7))
####    plt.title('-Spherical- Crystallite size dist for {200}')
####    plt.xlabel('Crystallite Size -Nm-')
####    plt.ylabel('Frequency')
####    plt.grid()
####    plt.plot(xforlognormal, yforlognormal,color = 'k')
####    plt.axvline(x=D_0, color ='r', linestyle= '--',label='median')
####    plt.legend()
####    plt.show()  
