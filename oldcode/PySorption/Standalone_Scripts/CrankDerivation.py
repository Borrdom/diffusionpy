import casadi as cs
import numpy as np
import matplotlib.pyplot as plt

t=cs.SX.sym("t")
x=cs.SX.sym("x")
D=1E-13
mp=0.001
d=0.01
B=np.pi/4*d**2
L=0.1
V=B*L
rhopinfty=mp/V

A=mp/(2*B*(np.pi*D)**(1/2))
C=A/t**(1/2)*cs.exp(-x**2/(4*D*t))

# dCdt=cs.gradient(C,t)
# d2Cdx2=cs.gradient(cs.gradient(C,x),x)
# null=dCdt-D*d2Cdx2

# null_fun=cs.Function("null_fun",[t,x],[null])
# print(null_fun(2,4))
tvec=np.linspace(0,1000,100)

C_fun=cs.Function("C_fun",[t,x],[C])
z_vec=np.linspace(0,L/2,100)
fig,ax=plt.subplots()
ax.set_xlabel("t[min]")
ax.set_ylabel("$/rho_p$[g/m^3]")
[ax.plot(tvec/60,1000*C_fun(tvec,z_vec[i]).full()) for i,val in enumerate(z_vec)]