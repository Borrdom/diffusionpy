
import numpy as np
import time
#from epcsaftpy import component, pcsaft
from numba import njit, config
#config.DISABLE_JIT = True
from PyCSAFT_nue import SAFTSAC,vpure,etaiter,THFaktor,ares

T=298.15
p=1E5
npoint=12000
mi=np.asarray([1.20469,4.2102])
sigi=np.asarray([2.797059952,3.0741])
ui=np.asarray([353.95,219.92])
epsAiBi=np.asarray([2425.67,1890.72])
kapi=np.asarray([0.04509,0.006692])
N=np.asarray([1.,1.])
xi=np.asarray([0.5,0.5])
x1=np.linspace(0,1,npoint)
x2=1-x1
xi=np.vstack((x1,x2))
vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)

kij=None
kijAB=None
start=time.time_ns()
lngammai=np.asarray([SAFTSAC(T,vpures,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB).flatten() for i,val in enumerate(xi[0,:])])
end=time.time_ns()
print((end-start)/1E9)
start=time.time_ns()
etamix=np.asarray([etaiter(p,T,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB) for i,val in enumerate(xi[0,:])])
results=[ares(T,val,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB) for i,val in enumerate(etamix)]
end=time.time_ns()
print((end-start)/1E9)
start=time.time_ns()
lnGammai=np.asarray([1+xi[:,i]*THFaktor(T,vpures,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB) for i,val in enumerate(xi[0,:])])
end=time.time_ns()
print((end-start)/1E9)
print("check")



# Water = component('Water', ms = 1.2047, sigma = 2.797 , eps = 353.94, 
#                kappaAB = 0.045099, eAB =2425.7, sites = [0, 1, 1], Mw = 18.015)

# Butanol = component('Butanol', ms = 4.2102, sigma = 3.0741 , eps = 219.92, 
#                kappaAB = 0.006692, eAB = 1890.72, sites = [0, 1, 1], Mw = 74.123)


# mix2=pcsaft(Water+Butanol)
# lnphip2=mix2.logfugef(xi[:,0],T,p,"L")[0]
# start=time.time_ns()
# res2=np.asarray([mix2.logfugef(xi[:,i],T,p,"L")[0] for i in range(npoint)])
# end=time.time_ns()
# print((end-start)/1E9)

