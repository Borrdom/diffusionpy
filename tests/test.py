import os
os.add_dll_directory("c:/Program Files (x86)/Intel/oneAPI/compiler/latest/windows/redist/intel64_win/compiler")
from numbalsoda import lsoda_sig, lsoda
from numba import njit, cfunc
import numpy as np
import matplotlib.pyplot as plt

@cfunc(lsoda_sig)
def rhs(t, u, du, p):
    du[0] = u[0]-u[0]*u[1]
    du[1] = u[0]*u[1]-u[1]*p[0]

funcptr = rhs.address # address to ODE function
u0 = np.array([5.,0.8]) # Initial conditions
data = np.array([1.0]) # data you want to pass to rhs (data == p in the rhs).
t_eval = np.linspace(0.0,50.0,1000) # times to evaluate solution

# integrate with lsoda method
usol, success = lsoda(funcptr, u0, t_eval, data = data)
plt.plot(t_eval,usol[:,0])
plt.plot(t_eval,usol[:,1])
plt.show()

# integrate with dop853 method

# usol = solution
# success = True/False