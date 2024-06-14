import numpy as np
from scipy.integrate import solve_ivp
from numba import njit,config,prange,jit
import time
from .diffusion import Diffusion_MS,wegstein
from scipy.interpolate import interp1d
from .PCSAFT import supersaturation

# config.DISABLE_JIT = True
# @njit(['Tuple((f8[:,::1], f8[:,::1]))(f8, f8[::1], f8[::1],  i8[::1], i8[::1],i8[::1], f8[::1], f8[::1],f8[::1],f8[::1], f8[::1], f8[::1], f8[::1], f8[::1],f8[::1],f8[::1],f8[::1],f8[:,::1],f8[:,::1])'],cache=True)
def CNT(t,alpha,A,B,n,T,lnS):
    """Calculates the crystallization kinetics based on the classical nucleation  coupled with a simple crystal growth model
    Args:
        t (array_like): time /s
        alpha (array_like): crystal fraction           /-
        A (float): crystallization rate constant m^2/s
        B (float): interfacial tension parameter             /N/m
        n (float): growth order                         /-
        T(float): Temperature                 /K
        lnS(array_like): log of supersaturation  /-
    Returns:
        ndarray:   
        recrystallization rate       /- \n
        growth rate    /- \n
    """
    N = T*np.exp(-B/lnS**2/T**3)
    G = T*(1-1/np.exp(lnS))
    # G=T*lnS
    k=A*N*G
    dalphadt=k*(t)**(n-1)*(1-alpha)
    dalphadt[alpha>1]=0
    dalphadt[alpha<0]=0
    return dalphadt


# config.DISABLE_JIT = True
def crystallization_mode(ode,A,B,n,T,saftpar,wv_fun=None):
    """alter the ode function in diffusionpy.Diffusion_MS, to also solve the crystallization
    Args:
        wvinit (array_like): vector of the mass fractions of the mobile components
        ode (array_like): ode fuinction which is modified by the function for the rest see CNT
    Returns:
        array_like: new modified ode function with the same format as the input ode function
    """
    # CNT1=njit(['f8[::1](f8, f8[::1], f8,  f8, f8 , i8, f8[::1])'],cache=True)(CNT)
    crystallizes=np.where(saftpar['deltaHSL']>0)[0][0]
    def crystallization_ode(*args):
        """solves the genralized Maxwell model for relaxation"""
        nTH,nz_1=args[-2].shape
        t=args[0]
        x=args[1]
        mobiles=args[4]
        immobiles=args[5]
        wi0=args[8]
        nc=len(wi0[:,0])
        argsmod=list(args)
        wv=np.zeros((nTH,nz_1))
        wi=np.zeros((nc,nz_1))
        for i in range(nTH): wv[i,:]=x[(nz_1)*(i):(nz_1)*(1+i)]
        alpha=x[(nz_1)*(nTH):(nz_1)*(nTH+1)]
        dl_la = np.fmax(np.fmin((wi0[crystallizes]-alpha*wi0[crystallizes])/(1-alpha*wi0[crystallizes]),1),0)
        # print(alpha)
        
    
        wi[mobiles,...]=wv
        if immobiles.shape[0]>0.:
            wi[crystallizes,...]=(1-np.sum(wv,axis=0))*dl_la
            wi[immobiles[crystallizes!=immobiles],...]=(1-np.sum(wv,axis=0))*(1-dl_la)
            
        else:
            wi[crystallizes,...]=dl_la
            wi[mobiles,...]=wi[mobiles,...]/np.sum(wi[mobiles,...],axis=0)

        if wv_fun is not None:
            wiB=args[-1].copy()
            wiB[:,mobiles]=wv_fun(dl_la[-1]) 
            wiB[:,immobiles[crystallizes!=immobiles]]=(1-np.sum(wiB[:,mobiles],axis=1))[:,None]*(1-dl_la[-1])
            wiB[:,crystallizes]=(1-np.sum(wiB[:,mobiles],axis=1))*dl_la[-1]
            argsmod[-1]=wiB
            wi[:,-1]=wiB[-1,:]
        
        lnS=np.asarray([supersaturation(T,val,**saftpar)[crystallizes] for val in wi.T])

        argsmod[1]=wv.flatten()

        dwvdt=ode(*argsmod)
        dalphadt=CNT(t,alpha,A,B,n,T,lnS)
        # fvec=np.hstack((dwvdt.flatten(),dalphadt.flatten(),drdt.flatten()))
        fvec=np.hstack((dwvdt.flatten(),dalphadt.flatten()))
        return fvec
    return crystallization_ode
