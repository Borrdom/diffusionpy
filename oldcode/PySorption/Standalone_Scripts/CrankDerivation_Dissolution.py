import casadi as cs
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
t=cs.SX.sym("t")
x=cs.SX.sym("x")
D=1E-10
mp=0.0005
d=0.01
B=np.pi/4*d**2
L=0.045
V=B*L #3.5 ml
rhow=997
rhop=1180
Vp=mp/rhop
rp=Vp/B#(Vp*3/4/np.pi)**(1/3)
mw=rhow*V
rhopinfty=mp/V

A=rhopinfty*L/2*(np.pi*D*t)**(-1/2)
C=A*cs.exp(-x**2/(4*D*t))
#C=rhop/2*(erf((-x+rp)/(2*(t*D)**(1/2)))+erf((x+rp)/(2*(t*D)**(1/2))))
C=sum([rhop/2*(erf((-x+rp+2*L*n)/(2*(t*D)**(1/2)))+erf((x+rp-2*L*n)/(2*(t*D)**(1/2)))) for n in range(-20,20)])

dCdt=cs.gradient(C,t)
d2Cdx2=cs.gradient(cs.gradient(C,x),x)
null=dCdt-D*d2Cdx2

null_fun=cs.Function("null_fun",[t,x],[null])
print(null_fun(2,4))
t0=0
tend=10000000
tvec=np.linspace(t0,tend,10000)

C_fun=cs.Function("C_fun",[t,x],[C])
z_vec=np.linspace(0,L,100)
fig,ax=plt.subplots()

ax.set_xlabel("t/min")
ax.set_ylabel(r'$\rho_p$/$kg/m^3$')
[ax.plot(tvec/60,C_fun(tvec,z_vec[i]).full()) for i,val in enumerate(z_vec)]
ax.plot([t0/60,tend/60],[rhopinfty,rhopinfty],'k--')