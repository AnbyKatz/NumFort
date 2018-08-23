#!/path/to/anaconda3/bin/python3.6

import matplotlib.pyplot as pypl
import numpy as np
import math as m
import pprint
from sympy import*
from matplotlib import rc

data = np.loadtxt("data.dat")
x = data[:,0]
y = data[:,1]
# z = data[:,2]
# w = data[:,3]

fig = pypl.figure()
axes = pypl.gca()

pypl.plot(x,y,label="plot 1")
pypl.xlabel("",fontsize=17)
pypl.ylabel("",fontsize=17)
pypl.title("" ,fontsize=23)
pypl.grid()

# pypl.plot(z,w,label="plot 2")
# pypl.legend()
# axes.set_xlim([8,16])
# axes.set_ylim([8,16])
# fig.savefig('graph.png')

pypl.show()
