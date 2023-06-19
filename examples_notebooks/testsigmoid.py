import numpy as np
import matplotlib.pyplot as plt
wi0=1E-4
wi8=0.1
a10=0.001
a18=0.25
dlna1_dlnw1=(np.log(a10)-np.log(a18))/(np.log(wi0)-np.log(wi8))
EJ=1E10
v2=1./997.
rho0=1180
R=8.31445
T=298.15
M2=18.015/1000.
RV=R*T/v2/M2
A=RV/EJ/v2
rhow0=wi0/(1-wi0)*rho0
rhow8=wi8/(1-wi8)*rho0
m=(a10-a18)/(rhow8-rhow0)
a_fun=lambda rho: a18+m*(rhow8-rho)
sig=lambda EJ : rhow0*(1+np.tanh(EJ/RV*v2))
EJvec=10**np.linspace(5,14,100)
plt.plot(np.log10(EJvec),sig(EJvec))
plt.show()