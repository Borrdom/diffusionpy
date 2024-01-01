import numpy as np
wvec=np.asarray([0,1,2])
gammavec=np.asarray([1,2,3])
import casadi as cs
gamma_fun=cs.interpolant('gamma_fun','linear',[wvec],gammavec)
import matplotlib.pyplot as plt
#from scipy.interpolate import InterpolatedUnivariateSpline
#gamma_fun=InterpolatedUnivariateSpline(wvec,gammavec)
plt.plot(wvec,gammavec,'kx')

wvec2=np.linspace(0,2,100)
gammavec2=gamma_fun(wvec2)
plt.plot(wvec2,gammavec2)
w=cs.MX.sym("w")
gamma=gamma_fun(w**2)