import numpy as np
from numba import njit
def averaging(a,axis=0):
    s1, e1 ,s2, e2 = 0,-1, 1, None 
    slc1 = [slice(None)] * a.ndim
    slc2 = [slice(None)] * a.ndim
    slc1[axis] = slice(s1, e1)
    slc2[axis] = slice(s2, e2)
    return (a[tuple(slc1)]+a[tuple(slc2)])/2.
b=np.ones((20,100)) 
c=averaging(b,axis=1)
print(c.shape)

@njit
def averaging2(a,axis=0):
    if axis==0: return (a[1:,:]+a[:-1,:])/2.
    if axis==1: return (a[:,1:]+a[:,:-1])/2.

d=averaging2(b,axis=1)
print(d.shape)

nc=3
nz=20
ji=np.zeros((nc,nz))
dji=np.diff(np.hstack((np.zeros((nc,1)),ji)))
print(dji.shape)


import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

x=np.linspace(0,1,11)

y=np.vstack((np.sin(x),np.cos(x),np.log1p(x)))

y=np.stack((y,y))

xi=np.linspace(0,1,101)

intf= interp1d(x,y,axis=2)

yi=intf(xi)

plt.ioff()
plt.plot(x,y.T,'x',
         xi,yi.T)
plt.show()
