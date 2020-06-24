import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt
from pymatgen.core.spectrum import Spectrum

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


plt.style.use(['seaborn-paper', 'presentation'])
plt.rc('font', family='serif')

material_list = ['asdep', '400', '710', '960', '1060']
tick_list     = [material + '$^\circ$C' for material in material_list]
tick_list[0]  = 'asdeposited'
two_theta = np.linspace(30,100,num=1000)
plt.figure(figsize = (14,10))
for c,material in enumerate(material_list):

    measurement = pd.read_excel('..//excelSheets//GADDS_{}.xlsx'.format(material),
                  sheet_name='RawSample')
    print(measurement)
    
    y = np.zeros(len(two_theta)) 

    for LP in range(0,6):
        y += pseudo_voigt(two_theta,measurement['peak_two_thetas'][LP],
                      measurement['fwhm'][LP],
                      measurement['fraction'][LP],
                      100)
    spect = Spectrum(two_theta,y)
    spect.normalize(mode="max", value=100)
    plt.plot(spect.x, spect.y + 110*c, c = 'k')
    plt.yticks([c*110 for c in range(len(tick_list))], tick_list,fontsize=20)
plt.xlabel('2\u03F4 $(^\circ)$')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='k', linestyle='--', alpha = 0.4)
plt.show()
