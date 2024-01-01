import casadi as cs
import numpy as np
from scipy.interpolate import CubicSpline
Psibin=np.asarray([0.044339926,0.158867985,0.276410999,0.42563113,0.618467599,1])
D=-np.asarray([1.25E+00,1.26E+00,1.32E+00,1.30E+00,1.28E+00,1.23E+00])
Psibin2=np.asarray([0.144339926,0.28867985,0.376410999,0.52563113,0.718467599,1])
D2=-(np.asarray([1.25E+00,1.26E+00,1.32E+00,1.30E+00,1.28E+00,1.23E+00])+0.1)

D4_fun=CubicSpline(Psibin,D,bc_type=((2, 1.0), (2, 1.0)))
D42_fun=CubicSpline(Psibin,D2,bc_type=((2, 1.0), (2, 1.0)))
xvec=np.linspace(0,1.5,100)
Dvec4=D4_fun(xvec)
Dvec42=D42_fun(xvec)

D_fun=cs.interpolant('LUT','linear',[xvec],Dvec4)
D2_fun=cs.interpolant('LUT','linear',[xvec],Dvec42)
bD_fun=cs.interpolant('LUT','bspline',[Psibin],D)
bD2_fun=cs.interpolant('LUT','bspline',[Psibin2],D2)

Dvec=D_fun(xvec).full().T[0]
Dvec2=D2_fun(xvec).full().T[0]
bDvec=bD_fun(xvec).full().T[0]
bDvec2=bD2_fun(xvec).full().T[0]

DASD=(0.5/Dvec+0.5/Dvec2)**-1
bDASD=(0.5/bDvec+0.5/bDvec2)**-1
import matplotlib.pyplot as plt
plt.plot(xvec,Dvec,'k-')
plt.plot(xvec,Dvec2,'r-')
plt.plot(xvec,DASD,'gx')
plt.plot(Psibin,D,'kx')
plt.plot(Psibin2,D2,'ko')
plt.plot(xvec,bDASD,'g--')
plt.plot(xvec,bDvec,'k--')
plt.plot(xvec,bDvec2,'r--')
plt.plot(xvec,Dvec4,'o-')
plt.show()
