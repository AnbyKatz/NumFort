#!/path/to/anaconda3/bin/python3.6

import matplotlib.pyplot as pypl
import numpy as np
import math as m
import pprint
from sympy import*
from matplotlib import rc

data = np.loadtxt("filename")
x = data[:,0]
y = data[:,1]

fig = pypl.figure()
axes = pypl.gca()

pypl.plot(x,y,label="1")
pypl.xlabel("",fontsize=17)
pypl.ylabel("",fontsize=17)
pypl.title("" ,fontsize=23)
pypl.grid()

# pypl.legend()
# axes.set_xlim([,])
# axes.set_ylim([,])
# fig.savefig('graph.png')

pypl.show()
