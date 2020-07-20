'''
    This module reads an synchrotron data and does warren-averbach analysis.
    Data is provided as Q not as 2theta.
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

indi = 401

def Q_to_h3(profile):
    q1 = profile[2] - 2*pi/profile[0]*0.5*10
    q2 = profile[2] + 2*pi/profile[0]*0.5*10
    x = np.linspace(q1,q2, num = 4001)
    y = pseudo_voigt(x,profile[2],
                       profile[1],
                       profile[4],100)
    h3_x = x*profile[0]/2/pi/10
    profileq = Spectrum(h3_x,y)
    return profileq

def get_fourier_coefs(profile_spec, plot = False):

    x = profile_spec.x
    h = profile_spec.y
    
    H = rfft(h)
    frequency = rfftfreq(len(x),round(x[1]-x[0],6))
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
    return freq[:indi],coefs[:indi]

##def WA(profile_spec1, profile_spec2 ):
##    freq1, coefs1 = get_fourier_coefs(profile_spec1, plot = False)    
##    freq2, coefs2 = get_fourier_coefs(profile_spec2, plot = False)
##    
##    coefslog1 = np.log(coefs1)
##    coefslog2 = np.log(coefs2)
#### This plot will plot the fourier coefficients of a profile
##    plt.figure(figsize=(10,8))
##    msZn     = []
##    msStrain = [0]
##    count = 0
##    for i,j in zip(coefslog1, coefslog2):
##        if count <6:
##            plt.plot([1, 4], [i, j], 'ko--')
##            count+=1
##        msZn.append((i-j)/3/2/(np.pi**2))
##    for n,Zn in enumerate(msZn[1:]):
##        msStrain.append(Zn/(n+1)/(n+1))
##        
##    msZn[0]   = 0.0 # gets negative value yet so close to zero
##    rmsZn     = [np.sqrt(i) for i in msZn]
##    rmsStrain = [np.sqrt(i) for i in msStrain]
##    plt.ylabel('$\ln(A_{n})$')
##    plt.xlabel('$l^2$')
##
##    plt.grid()
##    plt.show()
##    print('First 10 avg. strain in' + '{}'.format('111') + 'direction:\n {}'.format(rmsStrain[:10])  )
##    print('First 10 avg. Zn in' + '{}'.format('111') + 'direction:\n {}'.format(rmsZn[:10])  )
##    print('*-*-*-*')
##
##    Asize = np.exp(np.array(msZn)*2*1*(np.pi**2) +np.log(coefs1))
##
#### This plot will plot the size coefficients of a profile    
##    plt.figure(figsize=(8,8))
####    plt.title('Size coefficients for {'+ '{}'.format(plane_dir)+ '}')
##    plt.xlabel('n', fontsize = 25)
##    plt.ylabel('As', fontsize = 25)
##    plt.scatter(range(indi),Asize,color = 'k')   
##    reg = LinearRegression().fit( np.arange(5).reshape(-1,1),Asize[:5])   
##    dummyx = np.linspace(0,-reg.intercept_/reg.coef_,num=50)
##    dummyy = reg.coef_[0]*dummyx + reg.intercept_
##    plt.plot(dummyx, dummyy, 'r--')
##    plt.scatter(-reg.intercept_/reg.coef_,0,s=35)
##
##    plt.grid()
##    plt.show()
##    
##    return rmsStrain
##if __name__ == '__main__':
##    a = 2.099
####    profile1 = [2.423, 0.495, 25.93, 29.426, 0.673]
####    profile2 = [1.21, 0.798, 51.942, 27.46, 0.969]
##    profile1 = [2.099, 0.6, 29.932, 76.787, 0.549]
##    profile2 = [2.099, 1.255, 59.97, 15.216, 0.886]
##    profile_spec1 = Q_to_h3(profile1)
##    profile_spec2 = Q_to_h3(profile2)
##    WA(profile_spec1, profile_spec2 )
##    
##if __name__ == '__main__':
##       
##    strain111, Larea111, Lvolume111 = WA(111)
##    a111 = a
##    strain200, Larea200, Lvolume200 = WA(200)
##    a200 = a
#### This part is about strain versus L plot
##    plt.figure(figsize = (14,6.5))
####    plt.title('Strain over length L')
##    plt.xlabel('L (Å)')
##    plt.ylabel('RMS Strain (%)')
##    plt.grid()
##    plt.axhline(y=wh_strain, color ='r', linestyle= '--',label='W-H')
##    LL    = [i*a200 for i in range(17)]
##    LL111 = [i*a111 for i in range(15)]
##    plt.xticks(LL111)
##    plt.plot(LL111[1:],[i*100 for i in strain111[1:15]],'o-',label='111')
##    plt.plot(LL[1:],[i*100 for i in strain200[1:17]],'o-',label='200')
##    plt.ylim([0.25,2.1])
##    plt.legend()
##    plt.show()
##    
#### This part is about crystallite size distribution
##    Larea     = Larea111
##    Lvolume   = Lvolume111
##    plane_dir = '111'
##    
##    D_0, sigma = solve_log_normal(Larea, Lvolume).x
##    print('D_0(median)    , ' + '{}'.format(plane_dir) + ' plane : {}'.format(round(D_0,4)))
##    print('Sigma(variance), ' + '{}'.format(plane_dir) + ' plane : {}'.format(round(sigma,4)))
##
##    print('*-*-*-*')
##    
##    print('<D>area   ' + '{}'.format(plane_dir) + ' direction: {}'.format(D_0*np.exp(5/2*np.log(sigma)*np.log(sigma))))
##    print('<D>volume ' + '{}'.format(plane_dir) + ' direction: {}'.format(D_0*np.exp(7/2*np.log(sigma)*np.log(sigma))))
##    print('*-*-*-*')
##
##    xforlognormal    = np.linspace(0.0001, 80, num = 1500)
##    yforlognormal = log_normal(D_0, sigma, xforlognormal)
##    
##    plt.figure(figsize=(14,7))
####    plt.title('-Spherical- Crystallite size dist for {}'.format(plane_dir))
##    plt.xlabel('Crystallite Size(Nm)')
##    plt.ylabel('Frequency')
##    plt.grid()
##    plt.plot(xforlognormal, yforlognormal,color = 'k')
##    plt.axvline(x=D_0, color ='r', linestyle= '--',label='Median')
##    plt.axvline(x=D_0*np.exp(5/2*np.log(sigma)*np.log(sigma)), color ='b', linestyle= '--',label='<D>area')
##    plt.axvline(x=D_0*np.exp(7/2*np.log(sigma)*np.log(sigma)), color ='g', linestyle= '--',label='<D>volume')
##    plt.legend()
##    plt.show()  
