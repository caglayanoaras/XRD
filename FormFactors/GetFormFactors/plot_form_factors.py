#plot atomic form factor figures

from os             import listdir
from os.path        import isfile, join
from scipy.optimize import curve_fit
import sys
import json
import numpy             as np
import matplotlib.pyplot as plt

plt.rc('font', family='serif')
plt.style.use(['seaborn-paper', 'presentation'])

with open('..//atomic_scattering_params_ict.json') as f:
    ATOMIC_SCATTERING_PARAMETERS = json.load(f)
with open('..//atomic_scattering_params_ions.json') as f:
    ATOMIC_SCATTERING_PARAMETERS_ION = json.load(f)
    
def form_factor(element = 'N', ion = True, s_array = np.arange(0,1,0.005) ):

    if ion:
        coeff = ATOMIC_SCATTERING_PARAMETERS_ION[element]
    else:
        coeff = ATOMIC_SCATTERING_PARAMETERS[element]

    coeffs = np.array(coeff)

    factor_array = []

    for s in s_array:
        s2 = s*s
        
        fs =  np.sum(coeffs[:,0]*np.exp(-coeffs[:,1]*s2))

        factor_array.append(fs)
    return  factor_array

#*-*-*-*-*
def get_dft_formfactors():
    
    wlist = np.linspace(0,2*np.pi,num=60) #Check jupyter notebooks for verification
    try:
        form_ti_dft = np.load('576-576-576_111.npy')
    except:
        print('Check the DFT calculation of Ti form factor!')
        sys.exit()

    try:
        form_al_dft = np.load('192-192-192_111.npy')
    except:
        print('Check the DFT calculation of Al form factor!')
        sys.exit()

    try:
        form_N100_dft = np.load('192-384-192_100.npy')
    except:
        print('Check the DFT calculation of N_100 form factor!')
        sys.exit()

    try:
        form_N111_dft = np.load('192-384-192_111.npy')
    except:
        print('Check the DFT calculation of N_111 form factor!')
        sys.exit()
    summ = 7+13+22
    summ_dft = np.abs(form_ti_dft).max() + np.abs(form_al_dft).max() + np.abs(form_N100_dft).max()
    
    al   = np.abs(form_al_dft) * summ / summ_dft
    ti   = np.abs(form_ti_dft) * summ / summ_dft
    n100 = np.abs(form_N100_dft) * summ / summ_dft
    n111 = np.abs(form_N111_dft) * summ / summ_dft

    return al, ti, n100, n111, wlist

#*-*-*-*-*
def plot_ti_formfactors(omega_cut):

    s_array  = np.arange(0,1,0.005)
    form_ti3 = form_factor(element = 'Ti', ion = True,  s_array = s_array )
    form_ti  = form_factor(element = 'Ti', ion = False, s_array = s_array )
    w_array  = s_array * 4 * np.pi

    plt.figure(figsize=(10,6))
    plt.plot(w_array, form_ti3, label = 'Ti+3 ICT')
    plt.plot(w_array, form_ti,  label = 'Ti ICT')
    plt.plot(w_array_dft, ti,   label = 'Ti DFT')

    plt.xlim([0,omega_cut])
    plt.ylim([0,23])
    plt.grid()
    plt.legend()
    plt.xlabel('4\u03C0sin(\u03B8) / \u03BB')
    plt.ylabel('Form factor')
    plt.show()
    
#*-*-*-*-*
def plot_al_formfactors(omega_cut):

    s_array  = np.arange(0,1,0.005)
    form_al3 = form_factor(element = 'Al', ion = True,  s_array = s_array )
    form_al  = form_factor(element = 'Al', ion = False, s_array = s_array )
    w_array  = s_array * 4 * np.pi

    plt.figure(figsize=(10,6))
    plt.plot(w_array, form_al3, label = 'Al+3 ICT')
    plt.plot(w_array, form_al,  label = 'Al ICT')
    plt.plot(w_array_dft, al,   label = 'Al DFT')

    plt.xlim([0,omega_cut])
    plt.ylim([0,14])
    plt.grid()
    plt.legend()
    plt.xlabel('4\u03C0sin(\u03B8) / \u03BB')
    plt.ylabel('Form factor')
    plt.show()    

#*-*-*-*-*
def plot_n_formfactors(omega_cut):

    s_array  = np.arange(0,1,0.005)
    form_n = form_factor(element = 'N', ion = False,  s_array = s_array )
    w_array  = s_array * 4 * np.pi

    plt.figure(figsize=(10,6))
    plt.plot(w_array, form_n,  label = 'N ICT')
    plt.plot(w_array_dft, n100,  label = 'N(100) DFT')
    plt.plot(w_array_dft, n111,  label = 'N(111) DFT')

    plt.xlim([0,omega_cut])
    plt.ylim([0,14])
    plt.grid()
    plt.legend()
    plt.xlabel('4\u03C0sin(\u03B8) / \u03BB')
    plt.ylabel('Form factor')
    plt.show()

def fit_func(w_array_dft, a,b,c,d,e,f,g,h,i,j):
    
    s = w_array_dft / 4/np.pi
    ssquare_list = np.multiply(s,s)
    
    res = a*np.exp(-b*ssquare_list) + c*np.exp(-d*ssquare_list) +\
          e*np.exp(-f*ssquare_list) + g*np.exp(-h*ssquare_list) +\
          i*np.exp(-j*ssquare_list) 
    return res

def fit_form_analytically(wlist, formfact ):
    try:
        popt, pcov   = curve_fit(fit_func, wlist,formfact, bounds=(-600,600), maxfev=10000)
    except:
        popt, pcov   = curve_fit(fit_func, wlist,formfact, bounds=(-1000,1000), maxfev=10000)

    plt.title('quick check')
    plt.plot(w_array_dft,fit_func(w_array_dft, *popt))
    plt.plot(w_array_dft, formfact)
    plt.show()
    
    res = []
    for i in range(0,10,2):
        sset = [popt[i], popt[i+1]]
        res.append(sset)
        
    return res
    
if __name__ == '__main__':

    al, ti, n100, n111, w_array_dft = get_dft_formfactors()
    

    plot_ti_formfactors(2*np.pi)
    plot_al_formfactors(2*np.pi)
    plot_n_formfactors(2*np.pi)

# ALREADY DONE
####    dft_dict = {'N' :fit_form_analytically(w_array_dft, n111 ),
####                'Ti':fit_form_analytically(w_array_dft, ti   ),
####                'Al':fit_form_analytically(w_array_dft, al   ),}
####
####    with open('..//atomic_scattering_params_dft.json', 'w') as fp:
####        json.dump(dft_dict, fp)
        

