#!/usr/bin/python3.5
import matplotlib.pyplot as pypl
import numpy as np
import math as m
from scipy.integrate import quad
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
from pylab import*
from sympy import*

data = np.loadtxt("output.dat")
with open ("titles.dat", "r") as myfile:
    titles = myfile.readlines()
    
A = np.shape(data)
counter = 1
for ii in range(0,int(A[1]/2)):
    x = data[:,ii*2]
    y = data[:,ii*2+1]

    if titles[3] != "nil\n":
        pypl.plot(x,y,label=titles[ii+2])
    else:
        pypl.plot(x,y)
        
    counter = counter + 1

if titles[3] != "nil\n":
    pypl.legend()
    
pypl.xlabel(titles[0],fontsize = 18)
pypl.ylabel(titles[1],fontsize = 18)
pypl.grid()
pypl.show()
