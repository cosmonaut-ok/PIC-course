#!/usr/bin/env python3
import matplotlib.pyplot as plt
from numpy import genfromtxt,arange
Ey=genfromtxt("field_evol_Ey.dat")
Ey=Ey[:,1]
Bz=genfromtxt("field_evol_Bz.dat")
Bz=Bz[:,1]
arg=arange(len(Ey))
plt.plot(arg,Ey,label="$E_y$")
plt.plot(arg,Bz,label="$B_z$")
plt.xlabel("x")
plt.legend()
plt.savefig("field_evol_simple.pdf")
plt.show()
