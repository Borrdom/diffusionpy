import numpy as np
from scipy.optimize import root
from scipy.interpolate import interp1d
try:
    from numba import njit
except ImportError:
    def njit(f=None, *args, **kwargs):
        def decorator(func):
            return func

        if callable(f):
            return f
        else:
            return decorator

def time_dep_surface(t,wi0,wi8,mobile,taui,lngi_t=None):
    """_summary_

    Args:
        t (array_like): time vector
        wi0 (array_like): Mass fractions at t=0               
        wi8 (array_like): Mass fraction at t=infinity         
        mobile(array_like): boolean vector indicating the mobility of a component
        taui (array_like): time constant of the surface concentration function
        lngi_fun (array_like, optional): activity coefficient function.

    Returns:
        array_like: mass fraction at the surface as a vector of time
    """    
    nc=len(wi0)
    nf=np.sum(mobile)
    nt=len(t)
    allflux=nc==nf
    mobiles=np.where(mobile)[0] if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] if not allflux else np.asarray([-1],dtype=np.int64)
    wv8=wi8[mobiles]
    wv0=wi0[mobiles]
    lam=1-np.exp(-t[:,None]/taui[None,:])
    wvt=wv0*(1-lam)+wv8*lam

    if lngi_t is not None:
        lnai8=lngi_t[-1,:] +  np.log(wi8)
        lnai0=lngi_t[0,:]  +  np.log(wi0)
        av0=np.exp(lnai0[mobiles])
        av8=np.exp(lnai8[mobiles])
        avt=av0*(1-lam)+av8*lam
        lnwvt=np.log(avt)-lngi_t[:,mobiles]
        wvt=np.exp(lnwvt)
    wit=np.zeros((nt,nc))
    wit[:,mobiles]=wvt
    wit[:,immobiles]=(1-np.sum(wvt,axis=1)[:,None])*wi0[None,immobiles]/np.sum(wi0[immobiles])
    return wit