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