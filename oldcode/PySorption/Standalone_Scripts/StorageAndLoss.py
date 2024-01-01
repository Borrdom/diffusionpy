import numpy as np
import cmath
from numpy.random import rand
import matplotlib.pyplot as plt
nArm=2
logomega=np.linspace(-3,3,300)
omega=10**logomega
tauJ=10**np.linspace(-12,5,nArm)
tauJ[0]=1
tauJ[1]=1E2
EJ=rand(nArm)*(1E6-1E4)+1E4
EJ[0]=1E9
EJ[1]=1E5
Eeq=0
etaJ=tauJ*EJ

#pK=tauJ
#qK=EJ

#Num=sum([0+val*(omega*1j)**k for k,val in enumerate(pK)])
#Den=sum([0+val*(omega*1j)**k for k,val in enumerate(qK)])

#E=Num/Den
E=sum([Eeq+(omega*1j)/(1+val*omega*1j)*etaJ[k] for k,val in enumerate(tauJ)])

storage=E.real
loss=E.imag
#tandelta=storage/loss
etastar=np.sqrt((storage/omega)**2+(loss/omega)**2)
fig,ax=plt.subplots()
ax.plot(logomega,storage,'r-')
ax.plot(logomega,loss,'b-')
ax.plot(logomega,etastar,'g-')
plt.show()