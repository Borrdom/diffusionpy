import numpy as np
from numba import njit
from scipy.optimize import root
from scipy.interpolate import interp1d


def time_dep_surface(t,wi0,wi8,mobile,taui,lngi_fun=None):
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