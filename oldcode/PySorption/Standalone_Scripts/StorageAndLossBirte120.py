
import numpy as np
import cmath
from numpy.random import rand
import matplotlib.pyplot as plt
from  PySorption.DiffusionsmodellPlots import Plotlist
nArm=1
logomega=np.linspace(-5,5,300)
omega=10**(logomega)
tauJ=0.04525


EJ=4E9
Eeq=0
etaJ=tauJ*EJ
#pK=tauJ
#qK=EJ

#Num=sum([0+val*(omega*1j)**k for k,val in enumerate(pK)])
#Den=sum([0+val*(omega*1j)**k for k,val in enumerate(qK)])

#E=Num/Den7
#Maxwell
E=Eeq+(omega*1j)/(1+tauJ*omega*1j)*etaJ
#Voigt
Estorage=(omega**2)/(1/tauJ+tauJ*omega**2)*etaJ
Estorage=(tauJ**2*omega**2)/(1+tauJ**2*omega**2)*EJ


storage=E.real
loss=E.imag
print(Estorage-storage)
#tandelta=storage/loss
fig,ax=plt.subplots()
ax.plot(logomega,storage,'r-')
ax.plot(logomega,loss,'b-')
plt.show()
t=10**(-logomega)
#Plotlist([t],[storage],[t],[loss],Origin=True,xla="t/s",yla="E/Pas",filename="StoragePlot")
#ax.plot(logomega,tandelta,'g-')