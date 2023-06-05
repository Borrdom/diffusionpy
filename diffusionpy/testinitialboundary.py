import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
wi0=1E-4
wi8=0.1
a10=0.001
a18=0.25
dlna1_dlnw1=(np.log(a10)-np.log(a18))/(np.log(wi0)-np.log(wi8))
EJ=1E10
v2=1./997.
rho0=1180.
R=8.31448
T=298.15
M2=18.015/1000.
RV=R*T/v2/M2
A=RV/EJ/v2
B=A/rho0*dlna1_dlnw1
C=np.exp(np.log(wi8)*B)
print(B)
print(C)
rhow0=wi0*rho0
rhow8=wi8*rho0
m=(a10-a18)/(rhow8-rhow0)
a_fun=lambda rho: a18+m*(rhow8-rho)

deltaa_fun=lambda rho: dlna1_dlnw1*(-np.log(wi8)+np.log(rho/(rho0)))


rhow_plus= lambda rho: rhow0+A*np.log(a18/a_fun(rho))-rho
rhow_plus= lambda rho: rhow0-A*deltaa_fun(rho)-rho

rhovec=np.linspace(rhow0,rhow8,100000)
plt.plot(rhovec,rhow_plus(rhovec))

#plt.plot(rhovec[np.argmin(np.abs(rhow_plus(rhovec)))],np.min(np.abs(rhow_plus(rhovec))),"kx")
sol=rhovec[np.argmin(np.abs(rhow_plus(rhovec)))]
sigma_plus=RV*deltaa_fun(sol)
sigma_plus=RV*(np.log(a18/a_fun(sol)))
wsol=sol/(sol+rho0)
wsol_lamW=B*lambertw(wi8/B)
print(wsol)
print(wsol_lamW)
print(wsol_lamW*rho0)
sol2=dlna1_dlnw1*A
print(sol)
print(sigma_plus)
wsol2=sol2/(sol2+rho0)
print(wsol2)
plt.show()