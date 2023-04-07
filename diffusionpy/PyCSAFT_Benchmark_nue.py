
import numpy as np
from scipy.optimize import fsolve,root
import time
from PyCSAFT_nue import Z_obj,lnphii,SAFTSAC
from epcsaftpy import component, pcsaft

class mixture:
    def __init__(self,mi,sigi,ui,epsAiBi,kapi,N):
        self.mi=mi
        self.sigi=sigi
        self.ui=ui
        self.epsAiBi=epsAiBi
        self.kapi=kapi
        self.N=N

    def lnphi_p(self,T,p,xi):
        mi=self.mi
        sigi=self.sigi
        ui=self.ui
        epsAiBi=self.epsAiBi
        kapi=self.kapi
        N=self.N
        def TP_obj(x):
            eta=x
            return Z_obj(p,T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N)
        eta0=0.6
        eta_n=root(TP_obj,eta0, method='lm')["x"]
        return lnphii(T,eta_n[0],xi,mi,sigi,ui,epsAiBi,kapi,N)

        
T=298.15
p=1E5
npoint=12000
#Water Butanol
XAi=np.asarray([0.5,0.5])[:,None]
mi=np.asarray([1.2047,4.2102])[:,None]
sigi=np.asarray([2.797,3.0741])[:,None]
ui=np.asarray([353.94,219.92])[:,None]
epsAiBi=np.asarray([2425.7,1890.72])[:,None]
kapi=np.asarray([0.045099,0.006692])[:,None]
N=np.asarray([1.,1.])[:,None]

mix1=mixture(mi,sigi,ui,epsAiBi,kapi,N)
x1=np.linspace(0,1,npoint)
x2=(1-x1)

xi=np.stack((x1,x2))
# lnphip1=mix1.lnphi_p(T,p,xi[:,0][:,None])
# start=time.time_ns()
# res1=np.asarray([mix1.lnphi_p(T,p,xi[:,i][:,None]).flatten() for i in range(npoint)])
# end=time.time_ns()
# print((end-start)/1E9)

vpure=np.asarray([1000./18.015/1000.,789/46.07/1000.])
SAFTSAC(T,vpure,xi[:,0][:,None],mi,sigi,ui,epsAiBi,kapi,N).flatten()
start=time.time_ns()
res1=np.asarray([SAFTSAC(T,vpure,xi[:,i][:,None],mi,sigi,ui,epsAiBi,kapi,N).flatten() for i in range(npoint)])
end=time.time_ns()
print((end-start)/1E9)

#0.06

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
# print(res1-res2)

