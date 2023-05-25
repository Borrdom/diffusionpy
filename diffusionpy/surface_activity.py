import numpy as np
from numba import njit
@njit
def time_dep_surface(t,wi0,wi8,mobiles,immobiles,taui,rho):
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    rhoiB=np.zeros_like(wi0)
    rhoiB[mobiles]=(wi8[mobiles]+(wi0[mobiles]-wi8[mobiles])*np.exp(-t/taui))*rho if taui[0]!=0. else wi8[mobiles]*rho
    rhoiB[immobiles]=(rho-np.sum(rhoiB[mobiles],axis=0))*wi0_immobiles if taui[0]!=0. else rho*wi8[immobiles]
    return rhoiB