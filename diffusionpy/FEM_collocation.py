from numba import njit
import numpy as np

@njit(['Tuple((i8,i8,f8[:,::1],f8[::1]))(i8)'],cache=True)
def collocation_matrices(nz_1):
    nP=4
    nE=(nz_1-1)//(nP-1)
    NS=[0.0, 0.2928932188134524, 0.7071067811865476, 1.0]
    S=[[-5.828427124746187, 8.242640687119275, -3.414213562373087, 0.9999999999999973], [-1.4142135623730971, -0.41421356237309154, 2.414213562373092, -0.5857864376269044], [0.5857864376269034, -2.4142135623730896, 0.414213562373088, 1.4142135623730976], [-0.9999999999999956, 3.414213562373087, -8.242640687119275, 5.828427124746187]]
    return nP,nE,np.asarray(S),np.asarray(NS)

@njit(['f8[::1](f8[::1],f8[:,::1],i8,i8)'],cache=True)
def collocation_forward(x,S,nE,nP):
    dx=np.zeros_like(x)
    for i in range(nE):
        P0=i*(nP-1)
        P8=(i+1)*(nP-1)+1
        xP=x[P0:P8]
        if i==0: dxP_=np.zeros_like(xP)
        dxP=S@xP
        if i>0: dxP[0]=dxP_[-1]
        dx[P0:P8]=dxP
        dxP_=dxP
    return dx
@njit(['f8[::1](f8[::1],f8[:,::1],i8,i8)'],cache=True)
def collocation_backward(x,S,nE,nP):
    dx=np.zeros_like(x)
    for i in range(nE-1,-1,-1):
        P0=i*(nP-1)
        P8=(i+1)*(nP-1)+1
        xP=x[P0:P8]
        if i==(nE-1): dxP_=np.zeros_like(xP)
        dxP=S@xP
        if i<(nE-1): dxP[-1]=dxP_[0]
        dx[P0:P8]=dxP
        dxP_=dxP
    return dx

@njit(['f8[::1](f8[::1],i8,b1)'],cache=True)
def collocation(x,nz_1,forward):
    nP,nE,S,NS=collocation_matrices(nz_1)
    if forward:
        return collocation_forward(x,S,nE,nP)
    else:
        return collocation_backward(x,S,nE,nP)

@njit(['Tuple((f8[::1],i8))(i8)'],cache=True)
def collocation_space(nz_1):
    nP,nE,S,NS=collocation_matrices(nz_1)
    z=np.asarray([0.])
    zE=np.linspace(0,1,nE+1)
    for i in range(nE):
        z0=zE[i]
        z8=zE[i+1]
        zP=NS*(z8-z0)+z0
        z=np.hstack((z,zP[1:]))
    return z,nE 
