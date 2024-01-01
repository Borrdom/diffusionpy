import numpy as np
import cmath
from numpy.random import rand
import matplotlib.pyplot as plt
nArm=2
logomega=np.linspace(-10,15,300)
omega=10**logomega
tauJ=10**np.linspace(-12,5,nArm)
tauJ[0]=6.94*10**(0)
tauJ[1]=0.011*10**(0)
EJ=rand(nArm)*(1E6-1E4)+1E4
EJ[0]=3.64E10/3
EJ[1]=5.05E11/3
Eeq=0
etaJ=tauJ*EJ
etaEq=0
#pK=tauJ
#qK=EJ
#viscosity is the same as the spring constant if one maxwell arm. Else it can be much bigger or smaller
#Num=sum([0+val*(omega*1j)**k for k,val in enumerate(pK)])
#Den=sum([0+val*(omega*1j)**k for k,val in enumerate(qK)])
etawater=0.01


#E=Num/Den
EP=sum([Eeq+(omega*1j)/(1+val*omega*1j)*etaJ[k]+etaEq*omega*1j  for k,val in enumerate(tauJ)])
EW=etawater*omega*1j
phiW=0.01
phiP=1-phiW
E=phiP**2*EP+phiW**2*EW+2*phiW*phiP*np.sqrt(EW*EP)
#E=EP**phiP*EW**phiW


storage=E.real
loss=E.imag
viscosity=np.sqrt((storage/omega)**2+(loss/omega)**2)
viscosityP=np.sqrt((EP.real/omega)**2+(EP.imag/omega)**2)
#tandelta=storage/loss
fig,ax=plt.subplots()
fig1,ax1=plt.subplots()
fig2,ax2=plt.subplots()
ax.plot(logomega,storage,'r-')
ax.plot(logomega,loss,'b-')
ax.plot(logomega,viscosity,'g-')

ax.plot(logomega,EP.real,'r--')
ax.plot(logomega,EP.imag,'b--')
ax.plot(logomega,viscosityP,'g--')

ax2.plot(logomega,np.log10(storage),'r-')
ax2.plot(logomega,np.log10(loss),'b-')
ax2.plot(logomega,np.log10(viscosity),'g-')

ax2.plot(logomega,np.log10(EP.real),'r--')
ax2.plot(logomega,np.log10(EP.imag),'b--')
ax2.plot(logomega,np.log10(viscosityP),'g--')


ax1.plot(logomega,np.log10(abs((storage-EP.real)/EP.real)),'r-')
ax1.plot(logomega,np.log10(abs((loss-EP.imag)/EP.imag)),'b-')
ax1.plot(logomega,np.log10(abs((viscosity-viscosityP)/viscosityP)),'g-')
plt.show()
#ax.plot(logomega,tandelta,'g-')