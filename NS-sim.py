#!/usr/bin/env python3
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x,y,z=genfromtxt("NS-sim.dat").T
ax.plot(x,y,z)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
