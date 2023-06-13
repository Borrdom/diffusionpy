import numpy as np
from .Stefan_Maxwell_segmental import Diffusion_MS_iter,D_Matrix,Diffusion_MS
from .PyCSAFT_nue import vpure
from .surface_activity import time_dep_surface2
import matplotlib.pyplot as plt

nc=2
L=0.016
wi0=np.asarray([0.9,0.1])
wi8=np.asarray([0.04,0.96])
Mi=np.asarray([18.015,65000])

mobile=np.asarray([True,False])
T=298.15
p=1E5
nt=300

# texp=np.asarray([ 0., 2790.,  4410.,  6570.,  9660., 11910., 15630., 17820., 22620., 27960.,32280., 36720., 41040., 97920.])
texp=np.asarray([ 0, 3600,   7200,  11580,  15000,  17760,  19800,  24600,  28200,31800,  34200,  37800,  41460,  45240,  48600, 106800])
# ww_exp=np.asarray([0.95,0.94981791, 0.94952408, 0.94898034, 0.94791589, 0.94813047,0.94697211, 0.94674376, 0.94614307, 0.94404476, 0.93323287,0.92318329, 0.90732405, 0.09724013])
ww_exp=np.asarray([0.9,0.89299418, 0.88460461, 0.88365329, 0.866351  , 0.86052165,0.85309549, 0.83375686, 0.80719971, 0.77584219, 0.75664025,0.71118538, 0.64861202, 0.56269132, 0.34671022, 0.03925336])
t=np.linspace(0.,texp[-1],nt)


kij=D_Matrix(np.asarray([-0.156]),nc)
par={"mi":np.asarray([1.20469,2420.99]),
"si": np.asarray([2.797059952,2.947]),
"ui" :np.asarray([353.95,205.27]),
"eAi" :np.asarray([2425.67,0.]),
"kAi":np.asarray([0.04509,0.02]),
"NAi":np.asarray([1.,653.]),
"Mi": Mi,
"kij":kij}

vpures=vpure(p,T,**par)
par["vpure"]=vpures




Dvec=np.asarray([5.6E-10])
taui=np.asarray([21711.02587])

witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=None)
wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,witB=witB)
plt.plot(t/60,wt[:,0],'g',label="Ideal and wiB(t)")
plt.xlabel('t/min')
plt.ylabel('$w_i$/-')


Dvec=np.asarray([5E-9])
# taui=np.asarray([35711.02587])
taui=np.asarray([21711.02587])

witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=None)
wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,witB=witB,T=T,par=par)
plt.plot(t/60,witB[:,0],"bo")
plt.plot(t/60,wt[:,0],'b',label="Non-ideal and wiB(t)")
plt.legend("")

plt.plot(texp/60,ww_exp,"gD")
ARDwiB=np.sum(np.abs(1-np.interp(texp,t,wt[:,0])/ww_exp))

# taui=np.asarray([85711.02587])
taui=np.asarray([62711.02587])
witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=par)
wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,witB=witB,T=T,par=par)
plt.plot(t/60,witB[:,0],"ro")
plt.plot(t/60,wt[:,0],'r',label="Non-ideal and aiB(t)")

ARDaiB=np.sum(np.abs(1-np.interp(texp,t,wt[:,0])/ww_exp))
print(ARDwiB)
print(ARDaiB)



plt.show()
