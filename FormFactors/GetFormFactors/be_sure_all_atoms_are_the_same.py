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

wlist = np.linspace(0,2*np.pi,num=60)
form_1 = np.load('576-576-576_100.npy')    
form_2 = np.load('192-192-192_100.npy')    
form_3 = np.load('192-192-384_100.npy')
form_4 = np.load('192-384-192_100.npy')
form_5 = np.load('384-192-576_100.npy')
form_6 = np.load('384-384-192_100.npy') 


plt.figure(figsize=(10,6))
plt.plot(wlist, np.abs(form_1), label = 'form_1')
plt.plot(wlist, np.abs(form_2), label = 'form_2')
plt.plot(wlist, np.abs(form_3), label = 'form_3')
plt.plot(wlist, np.abs(form_4), label = 'form_4')
plt.plot(wlist, np.abs(form_5), label = 'form_5')
plt.plot(wlist, np.abs(form_6), label = 'form_6')
plt.grid()
plt.legend()
plt.xlabel('4\u03C0sin(\u03B8) / \u03BB')
plt.ylabel('Form factor')
plt.show()
