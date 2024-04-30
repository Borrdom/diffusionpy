import matplotlib.pyplot as plt
import numpy as np
t1=np.linspace(0,2*np.pi,11)
t2=np.linspace(0,2*np.pi,101)
fx=lambda t: 16*np.sin(t)**3
fy=lambda t: 13*np.cos(t)-5*np.cos(2*t)-2*np.cos(3*t)-np.cos(4*t)

x1=fx(t1)
y1=fy(t1)
x2=fx(t2)
y2=fy(t2)
fig,ax=plt.subplots()
ax.plot(x2,y2,'k-')
ax.plot(x1,y1,'ko')
ax.plot(x2/1.3,y2/1.3,'C6-')
ax.plot(x1/1.3,y1/1.3,'C6s')
ax.plot(x2/1.7,y2/1.7,'C3-')
ax.plot(x1/1.7,y1/1.7,'C3^')
ax.plot(x2/2.2,y2/2.2,'C0-')
ax.plot(x1/2.2,y1/2.2,'C0*')
ax.set_xlabel(r'$m_{Aminoacid}/(gmol^{-1})$')
ax.set_ylabel(r'$\alpha/(Jmol^{-1})$')
ax.set_xticks(np.linspace(-14,14,5))
ax.set_yticks(np.linspace(-14,14,5))
plt.show()


from numba import njit
from diffusionpy.PCSAFT import dlnai_dlnxi,vpure
from diffusionpy.diffusion import D_Matrix

T=298.15
p=1E5
nc=3 # number of components
L=1E-4 # estimated thickness of the ASD
wi0=np.asarray([0.49999,0.49999,0.00002]) # Weight fractions of the components API, polymer, water at t=0
wi8=np.asarray([0.00001,0.00001,0.99998]) # Weight fractions of the components API, polymer, water at t=8.T
Mi= np.asarray([357.79,65000,18.15])
kij=D_Matrix(np.asarray([-0.0621,-0.025,-0.156]),nc)
par={"mi":np.asarray([14.283,2420.99,1.2046 ]),
"si": np.asarray([3.535,2.947, 2.797059952]),
"ui" :np.asarray([262.79,205.27,353.95 ]),
"eAi" :np.asarray([886.4,0.,2425.67 ]),
"kAi":np.asarray([ 0.02,0.02,0.04509 ]),
"NAi":np.asarray([3.,653., 1.]),
"Mi": Mi,
"kij":kij,
"kijA":np.asarray([[0.]])}
vpures=vpure(p,T,**par)
par["vpure"]=vpures
wi=np.asarray([0.49999,0.49999,0.00002])
mobile=np.asarray([True,True,False])
# @njit
def massbalancecorrection(T,wi,wi0,mobile,saftpar):
    Mi=saftpar["Mi"]
    ri=Mi/np.min(Mi)
    dlnai_dlnwi=dlnai_dlnxi(T,np.ascontiguousarray(wi),**saftpar)
    Gammaij=dlnai_dlnwi/wi[None,:]*wi[:,None]/ri[:,None]
    if not np.all(mobile):
        mobiles=np.where(mobile)[0] 
        immobiles=np.where(~mobile)[0]
        wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
        correction=np.sum(Gammaij[mobiles,:][:,immobiles]*wi0_immobiles[None,:],axis=-1)  
        THFaktor=(Gammaij[mobiles,:][:,mobiles]-correction[:,None]) #*wi[...,None,mobiles]
    else:
        THFaktor=Gammaij
    return THFaktor

massbalancecorrection(T,wi,wi0,mobile,par)




from numba import njit

@njit
def foo(a,b,c):
    return a+b,c

par={'c':1,'b':2,'a':3}

foo(**par)