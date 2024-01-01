from Pylomer.PyGT import Tggt
import numpy as np
import matplotlib.pyplot as plt
w1=np.linspace(0,1,100)
w2=1-w1
w3=np.zeros_like(w1)
rhoi=np.asarray([1320,997,1180])
Tgi=np.asarray([317.6,136,420])
wi=np.asarray([w1,w2,w3])


Tgw2=Tggt(wi.T,rhoi,Tgi,q=0)


delTg=np.subtract.outer(Tgi,Tgw2)

mxb=(1-wi)/delTg
m=np.gradient(mxb,axis=1)/np.gradient(wi,axis=1)
ms=np.nanmean(m,axis=1)
Lambda1=0.062
Lambda2=0.039
Xhi1=Lambda1*ms
Xhi2=Lambda2*ms

plt.plot(wi[0,:],mxb[0,:])
plt.plot(wi[1,:],mxb[1,:])
plt.show()
