import numpy as np
from numba import njit
@njit
def time_dep_surface(t,wi0,wi8,mobiles,immobiles,taui,rho):
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    rhoiB=np.zeros_like(wi0)
    rhoiB[mobiles]=(wi8[mobiles]+(wi0[mobiles]-wi8[mobiles])*np.exp(-t/taui))*rho if taui[0]!=0. else wi8[mobiles]*rho
    rhoiB[immobiles]=(rho-np.sum(rhoiB[mobiles],axis=0))*wi0_immobiles if taui[0]!=0. else rho*wi8[immobiles]
    return rhoiB

def gassided(t,rhov,wv0,wv8,taui,THFaktor_fun,tint):
    THFaktor_ave=np.mean(THFaktor_fun(tint),axis=0)
    #THFaktor_ave=np.eye(2)
    #THFaktor_80=(THFaktor_fun(tint[0])+THFaktor_fun(tint[-1]))/2
    nTH=len(taui)
    nz_1=len(rhov)//nTH
    rhov=np.ascontiguousarray(np.reshape(rhov,(nTH,nz_1))) 
    lnai8_lnai0=THFaktor_ave@(np.log(wv8)-np.log(wv0))
    #lnai8_lnai0=THFaktor_80@(np.log(wv8)-np.log(wv0))
    #lnwidt=np.linalg.solve(THFaktor_fun(t),-1/taui*np.exp(-t/taui)*lnai8_lnai0)

    lnwidt=np.linalg.solve(THFaktor_ave,-1/taui*np.exp(-t/taui)*lnai8_lnai0)
    lnwidt=-1/taui*np.exp(-t/taui)*(np.log(wv8)-np.log(wv0))
    drhovBdt= -rhov[:,-1]*lnwidt
    #drhovBdt=-(wv0-wv8)*np.exp(-t/taui)*1200*1/taui
    return drhovBdt