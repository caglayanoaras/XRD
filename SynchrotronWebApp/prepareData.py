import numpy as np
import matplotlib.pyplot as plt

from os      import listdir
from os.path import isfile, join

folder_path = '31_447_04_TiAlN'
files = [f for f in listdir(folder_path) if isfile(join(folder_path, f))]

# Get files ends with .subtr
measurement_file = [i for i in files if i[-6:] == '.subtr']

temperature_list           = []
measurement_list_q         = []
measurement_list_intensity = []

for measurement in measurement_file:
    print(measurement)
    
    f = open(folder_path + '/' + measurement, 'r')
    
    number_of_comment_line = 3
    for i in range(number_of_comment_line):
        comment = f.readline()
        if '#T=' in comment:
            title = comment
            
    data = np.genfromtxt(f, delimiter=' ', comments='#')
##    plt.title(title)
##    plt.plot(data[:,0],data[:,1])
##    plt.show()
    temperature_list.append(float(title[3:8]))
    measurement_list_q.append(data[:,0])
    measurement_list_intensity.append(data[:,1])

    f.close()

import pandas as pd
df_dict = {'temperature': temperature_list,
           'q': measurement_list_q,
           'intensity': measurement_list_intensity}

df = pd.DataFrame(df_dict)
df.to_pickle(folder_path + '_df.pickle' )
