import numpy as np
from .Stefan_Maxwell_segmental import Diffusion_MS_iter,D_Matrix,Diffusion_MS
from .PyCSAFT_nue import vpure
from .surface_activity import time_dep_surface2
import matplotlib.pyplot as plt

nc=3
L=0.001
wi0=np.asarray([0.01,0.495,0.495])
wi8=np.asarray([0.9,0.05,0.05])
Mi=np.asarray([18.015,65000,357.57])
T=298.15
p=1E5

kij=D_Matrix(np.asarray([-0.156,-0.025,-0.0621]),nc)
par={"mi":np.asarray([1.20469,2420.99, 14.283]),
"si": np.asarray([2.797059952,2.947, 3.535]),
"ui" :np.asarray([353.95,205.27, 262.79]),
"eAi" :np.asarray([2425.67,0., 886.4]),
"kAi":np.asarray([0.04509,0.02, 0.02]),
"NAi":np.asarray([1.,653., 3.]),
"Mi": Mi,
"kij":kij}

vpures=vpure(p,T,**par)

par["vpure"]=vpures
Dvec=np.asarray([1E-11,2.3E-11,1.7E-15])
nt=300
t=np.linspace(0,300,nt)*60
mobile=np.asarray([True,False,False])
wtid=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True)
# wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,T=T,par=par)
plt.plot(t/60,wtid[:,0])
plt.plot(t/60,wtid[:,1])
plt.plot(t/60,wtid[:,2])
# plt.plot(t/60,wt[:,0])
plt.show()