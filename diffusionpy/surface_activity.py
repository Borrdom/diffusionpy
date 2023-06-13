import numpy as np
from numba import njit
from .PyCSAFT_nue import lngi
from scipy.optimize import root
from scipy.interpolate import interp1d
@njit
def time_dep_surface(t,wi0,wi8,mobiles,immobiles,taui,rho):
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    rhoiB=np.zeros_like(wi0)
    rhoiB[mobiles]=(wi8[mobiles]+(wi0[mobiles]-wi8[mobiles])*np.exp(-t/taui))*rho if taui[0]!=0. else wi8[mobiles]*rho
    rhoiB[immobiles]=(rho-np.sum(rhoiB[mobiles],axis=0))*wi0_immobiles if taui[0]!=0. else rho*wi8[immobiles]
    return rhoiB

def gassided(t,rhov,wv0,wv8,taui,THFaktor_fun,tint):
    THFaktor_ave=np.mean(THFaktor_fun(tint),axis=0)
    THFaktor_80=np.mean(THFaktor_fun(tint[np.asarray([0,-1])]),axis=0)
    nTH=len(taui)
    nz_1=len(rhov)//nTH
    rhov=np.ascontiguousarray(np.reshape(rhov,(nTH,nz_1))) 
    lnai8_lnai0=THFaktor_ave@(np.log(wv8)-np.log(wv0))
    ai8_ai0=np.exp(lnai8_lnai0)
    lam=1-np.exp(-t/taui)
    dlamdt=-1/taui*np.exp(-t/taui)
    ait_ai0=(lam*(ai8_ai0-1)+1)
    lnwidt=np.linalg.solve(THFaktor_fun(t),dlamdt*(ai8_ai0-1)/ait_ai0)
    # lnwidt=np.linalg.solve(THFaktor_80,dlamdt*(ai8_ai0-1)/ait_ai0)
    drhovBdt= -rhov[:,-1]*lnwidt 
    #drhovBdt=-(wv0-wv8)*np.exp(-t/taui)*1/taui*1200#np.sum(rhov[:,-1],axis=0)
    return drhovBdt

def time_dep_surface2(t,wi0,wi8,mobile,taui,par=None):
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

    if par is not None:
        T=298.15
        lnai8=lngi(T,wi8,**par)+np.log(wi8)
        lnai0=lngi(T,wi0,**par)+np.log(wi0)
        
        av0=np.exp(lnai0[mobiles])
        av8=np.exp(lnai8[mobiles])
        avt=av0*(1-lam)+av8*lam
        def obj_f(wv,av):
            wi=np.zeros_like(wi0)
            wi[mobiles]=wv
            wi[immobiles]=(1-np.sum(wv))*wi0[immobiles]/np.sum(wi0[immobiles])
            return av/np.exp(lngi(T,wi,**par))[mobiles]-wv
        wvt=np.asarray([root(obj_f,x0=wvt[i,:],args=(val))["x"] for i,val in enumerate(avt)])
    wit=np.zeros((nt,nc))
    wit[:,mobiles]=wvt
    wit[:,immobiles]=(1-np.sum(wvt,axis=1)[:,None])*wi0[None,immobiles]/np.sum(wi0[immobiles])
    return wit