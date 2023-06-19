import numpy as np
from diffusionpy import lngi,vpure,dlnai_dlnxi,D_Matrix
import matplotlib.pyplot as plt
nc=3
Mi=np.asarray([18.015,25700,357.79])
T=298.15
p=1E5
npoint=100
mi=np.asarray([1.20469,1045.6,14.283])
sigi=np.asarray([2.797059952,2.71,3.535])
ui=np.asarray([353.95,205.599,262.791])
epsAiBi=np.asarray([2425.67,0.,886.4])
kapi=np.asarray([0.04509,0.02,0.02])
N=np.asarray([1.,231.,3.])
vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)
kij=D_Matrix(np.asarray([-0.148,0.,0.]),nc)
ww=np.linspace(0,1,npoint)
DL=0.5
wa=DL*(1-ww)
wp=(1-DL)*(1-ww)
wt=np.vstack((ww,wp,wa)).T
xt=np.asarray([wt[i,:]/Mi/np.sum((wt[i,:]/Mi)) for i in range(npoint)])
xw=xt[:,0]
xp=xt[:,1]
xa=xt[:,2]


Gammaiw=np.asarray([dlnai_dlnxi(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij=kij) for i in range(npoint)]).T
lngammaiw=np.asarray([lngi(T,vpures,np.ascontiguousarray(wt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,Mi,kij=kij)[1] for i in range(npoint)])
Gammaiw_=np.gradient(lngammaiw+np.log(ww),np.log(ww))

Gammaiwkorr=Gammaiw[0,0,:]-Gammaiw[0,1,:]*(1-DL)-Gammaiw[0,2,:]*DL

Gammaix=np.asarray([dlnai_dlnxi(T,vpures,np.ascontiguousarray(xt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,kij=kij) for i in range(npoint)]).T
lngammaix=np.asarray([lngi(T,vpures,np.ascontiguousarray(xt[i,:]),mi,sigi,ui,epsAiBi,kapi,N,kij=kij)[1] for i in range(npoint)])

Gammaix_=np.gradient(lngammaix+np.log(xw),np.log(xw))

#plt.plot(Gammaiwkorrr,,'kx')

#np.eye(3)-wt[1,:]
A1vec=[]
A2vec=[]

for i in range(npoint):
    A1=Gammaix[:,:,i]@(np.eye(nc)-xt[i,:]).T
    A1=(np.eye(nc)-xt[i,:])@Gammaix[:,:,i].T
    A2=Gammaiw[:,:,i]@(np.eye(nc)-wt[i,:]).T
    A2=((np.eye(nc)-wt[i,:])@Gammaiw[:,:,i].T)*np.outer(Mi,Mi**-1)
    A3=(np.eye(nc)-np.outer(xt[i,:]**-1,wt[i,:]))@Gammaiw[:,:,i].T*np.outer(Mi,Mi**-1)
    #A3vec.append(np.sum(Gammaix[:,:,i].T-A3,0))
    #A3=np.sum(Gammaix[:,:,i]@(np.zeros((nc,nc))-xt[i,:]).T)
    #A4=np.sum(Gammaiw[:,:,i]@(np.zeros((nc,nc))-wt[i,:]).T)

    #A1=(Gammaix[:,:,i]@np.eye(nc))@(np.eye(nc)-xt[i,:]).T
    #A2=(Gammaiw[:,:,i]@np.eye(nc))@(np.eye(nc)-wt[i,:]).T

    # A1=Gammaix[:,:,i]@(np.zeros((nc,nc))-xt[i,:]).T
    # A2=Gammaiw[:,:,i]@(np.zeros((nc,nc))-wt[i,:]).T
    # A1=Gammaix[:,:,i]@(np.ones((nc,nc))-xt[i,:]).T
    # A2=Gammaiw[:,:,i]@(np.ones((nc,nc))-wt[i,:]).T
    A1vec.append(A1)
    A2vec.append(A2)
    #plt.plot(wt[i,:],np.sum(Gammaix[:,:,i].T-A3,1),'kx')
    plt.plot(A1,A2,'kx')



#plt.plot(Gammaix_/Gammaiw_,(1-xw)/(1-ww)) #its true
#plt.plot(Gammaix[0,0,:]/Gammaiw[0,0,:],(1-xw)/(1-ww)*)

#plt.plot(Gammaiw[0,0,:],(Gammaix[0,0,:]-(Gammaix[0,0,:]+Gammaix[1,0,:]+Gammaix[2,0,:])*xw)/(1-ww),'kx')
#Mmean=np.sum(xt*Mi,axis=1)
#plt.plot(Gammaiw[0,0,:]*(1-ww),Gammaix[0,0,:],'kx')

#plt.plot(Gammaiw[0,0,:]-ww*Gammaiw[0,0,:]-wp*Gammaiw[0,1,:]-wa*Gammaiw[0,2,:],Gammaix[0,0,:]-xw*Gammaix[0,0,:]-xp*Gammaix[0,1,:]-xa*Gammaix[0,2,:],'kx')
plt.show()