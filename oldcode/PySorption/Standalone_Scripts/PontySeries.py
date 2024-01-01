import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
plt.close('all')
nArm=3
EJ=np.linspace(1E7,1E4,nArm)
#EJ=rand(nArm)*(1E6-1E3)+1E4
#EJ[3]=1E3
tauJ=10**np.linspace(-10,5,nArm)
tauJ[-1]=10**-0.5
etaJ=tauJ*EJ
logtvec=np.linspace(-12,5,300)
tvec=10**logtvec
F=lambda t: sum([val1*np.exp(-t/val2) for i,(val1,val2) in enumerate(zip(EJ,tauJ))])
Fvec=F(tvec)
fig,ax=plt.subplots()
fig1,ax1=plt.subplots()
ax.plot(np.log10(tvec),np.log10(Fvec),"k-")
[ax.plot([np.log10(tauJ[i]),np.log10(tauJ[i])],[np.log10(1E4),np.log10(1E7)],"k-") for i,val in enumerate(tauJ)]
ax1.plot(EJ)
