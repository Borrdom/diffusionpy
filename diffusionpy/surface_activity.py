import numpy as np
from numba import njit
from scipy.optimize import root
from scipy.interpolate import interp1d


def time_dep_surface(t,wi0,wi8,mobile,taui,lngi_fun=None):
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

    if lngi_fun is not None:
        lnai8=lngi_fun(wi8)+np.log(wi8)
        lnai0=lngi_fun(wi0)+np.log(wi0)
        
        av0=np.exp(lnai0[mobiles])
        av8=np.exp(lnai8[mobiles])
        avt=av0*(1-lam)+av8*lam
        def obj_f(wv,av):
            wi=np.zeros_like(wi0)
            wi[mobiles]=wv
            wi[immobiles]=(1-np.sum(wv))*wi0[immobiles]/np.sum(wi0[immobiles])
            return av/np.exp(lngi_fun(wi))[mobiles]-wv
        wvt=np.asarray([root(obj_f,x0=wvt[i,:],args=(val))["x"] for i,val in enumerate(avt)])
    wit=np.zeros((nt,nc))
    wit[:,mobiles]=wvt
    wit[:,immobiles]=(1-np.sum(wvt,axis=1)[:,None])*wi0[None,immobiles]/np.sum(wi0[immobiles])
    return wit