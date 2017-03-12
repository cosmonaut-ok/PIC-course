#!/usr/bin/env python3
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x,y,z=genfromtxt("pmove.dat").T
ax.plot(x,y,z)
plt.show()
