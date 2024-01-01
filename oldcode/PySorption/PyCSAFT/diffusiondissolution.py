import numpy as np
import casadi as cs
import matplotlib.pyplot as plt
nz=30
c=cs.SX.sym("c",nz)
t=cs.SX.sym("t")
L=10E-6
z=np.linspace(0,L,nz)
dz=L/nz
D=1E-13

Jbound=0
J=D*(c[1:]-c[:-1])/dz
J=cs.vertcat(Jbound,J,Jbound)
dJdz=(J[1:]-J[:-1])/dz
dcdt=dJdz#
c0=np.zeros(nz)
c0[0:nz//2]=10
c0[nz//2::]=0
def SolveODEs(x,f,x0,opts,Time):
    dae={"x":x,"t":Time,"ode":f}#,"print_stats":True,"abstol":1E-12}#,max_num_steps":1000"collocation_scheme":"legendre","interpolation_order":2,"number_of_finite_elements":100,"rootfinder_options":{"error_on_fail": False}
    fint=cs.integrator("rho2int","cvodes",dae,opts)
    F=fint(x0=x0)
    return F["xf"]
tspan=np.linspace(0,1000,100)
opts= {"grid":tspan,"max_num_steps":10000,"regularity_check":True,"output_t0":True}
csol=SolveODEs(c,dcdt,c0,opts,t)


[plt.plot(tspan,csol[i,:].full().flatten()) for i in range(nz)]
plt.show()

