#!/path/to/anaconda3/bin/pythonVersion

import matplotlib.pyplot as pypl
import numpy as np
import math as m
from matplotlib import rc

data = np.loadtxt("output.dat")
with open ("titles.dat", "r") as myfile:
    titles = myfile.readlines()
    
A = np.shape(data)
for ii in range(0,int(A[1]/2)):
    x = data[:,ii*2]
    y = data[:,ii*2+1]
    pypl.plot(x,y,label=titles[ii+3])

if int(A[1])/2 >= 2:
    pypl.legend()
    
pypl.xlabel(titles[1],fontsize = 18)
pypl.ylabel(titles[2],fontsize = 18)
pypl.title(titles[0])
pypl.grid()
pypl.show()
