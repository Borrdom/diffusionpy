
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
import xloil as xlo



@njit(cache=True)
def rayleigh(t,xi,gint,tint,p0i):
    gi=np.ones_like(xi)
    nc=len(xi)
    for i in range(nc):
        gi[i]=np.interp(t,tint,gint[:,i])
    yi=np.zeros_like(xi)
    yi=(xi*gi*p0i)/np.sum(xi*gi*p0i)
    dxidt=xi-yi
    for i in range(nc):
        if xi[i]<1E-20: dxidt[i]=0 
    return dxidt

if __name__=="__main__":
    from diffusionpy import lngi,vpure,D_Matrix
    T=303.15
    p=1E5
    nc=4 # number of components
    wi0=np.asarray([0.100001667,0.799998333,0.03,0.07])
    wi8=np.asarray([0.005951831,0.158725973,0.250596659,0.584725537])
    kijvec=np.asarray([-0.045000,-0.022000,-0.128000,-0.001000,0.002673,-0.062100])
    kijHBvec=np.asarray([0,0,0,0,-0.455211419,0])
    Mi=np.asarray([18.015,46.069,357.79,65000.])
    mi=np.asarray([1.2047,2.3827,14.283,2420.99])
    sigi=np.asarray([2.79533,3.1771,3.535,2.947])
    ui=np.asarray([353.95,198.24,262.791,205.27])
    epsAiBi=np.asarray([2425.7,2653.4,886.4,0.])
    kAiBi=np.asarray([0.045099,0.032384,0.02,0.02])
    Na=np.asarray([1.,1.,3.,653.])
    vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kAiBi,Na)
    kij=D_Matrix(kijvec,nc)
    kijHB=D_Matrix(kijHBvec,nc)
    par={"mi":mi,
    "si": sigi,
    "ui" :ui,
    "eAi" :epsAiBi,
    "kAi":kAiBi,
    "NAi":Na,
    "Mi": np.ones(nc),
    "kij":kij,
    "kijA":kijHB,
    "vpure":vpures}

    gi_fun=lambda xi: np.exp(lngi(T,xi,**par))
    nt=100
    tint=np.linspace(0,10,nt)
    xi0=wi0/Mi/np.sum(wi0/Mi)
    p0i=np.asarray([42470,104000,0,0])

    gint=np.asarray([np.ones_like(xi0)]*nt)
    for i in range(10):
        xsol=odeint(rayleigh,xi0,tint,args=(gint,tint,p0i),tfirst=True)
        gint=np.asarray([gi_fun(val) for val in xsol])
        wsol=xsol*Mi/np.sum(xsol*Mi,axis=1)[:,None]
        plt.plot(tint,wsol)
        # plt.plot(tint,wsol[:,1])
    import pandas as pd
    pd.DataFrame(wsol).to_clipboard(header=None)
    plt.show()
@xlo.func
def rayleigh_xloil(t:xlo.Array(float,dims=1),wi0:xlo.Array(float,dims=1),gint:xlo.Array(float,dims=2),p0i:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1)):
    xi0=wi0/Mi/np.sum(wi0/Mi)   
    xsol=odeint(rayleigh,xi0.copy(),tint.copy(),args=(gint.copy(),tint.copy(),p0i.copy()),tfirst=True)
    wsol=xsol*Mi/np.sum(xsol*Mi,axis=1)[:,None]
    return wsol
