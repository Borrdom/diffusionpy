import numpy as np
import cmath
from numpy.random import rand
import matplotlib.pyplot as plt
nArm=1
logt=np.linspace(-5,5,300)
t=10**logt
epsilon0=1
epsilon=t*epsilon0
tauJ=10**0

EJ=1000

etaJ=tauJ*EJ
#pK=tauJ
#qK=EJ

#Num=sum([0+val*(omega*1j)**k for k,val in enumerate(pK)])
#Den=sum([0+val*(omega*1j)**k for k,val in enumerate(qK)])

#E=Num/Den
#E=sum([Eeq+(omega*1j)/(1+val*omega*1j)*etaJ[k] for k,val in enumerate(tauJ)])
sigma1=tauJ*epsilon0*EJ*(1-np.exp(-t/tauJ))
sigma2=epsilon0*EJ*(np.exp(-t/tauJ))
E1=sigma1/epsilon
E2=sigma2/epsilon0
#storage=E.real
#loss=E.imag
#tandelta=storage/loss
fig,ax=plt.subplots()
ax.plot(logt,sigma1,'r-')
ax.plot(logt,sigma2,'b-')
#ax.plot(logomega,loss,'b-')
#ax.plot(logomega,tandelta,'g-')
plt.show()