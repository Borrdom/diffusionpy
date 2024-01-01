import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

t=np.asarray(t)
m=np.asarray(m)
tvec=np.linspace(t[0]**(1/2),t[-1]**(1/2),300)**2
mvec=InterpolatedUnivariateSpline(t,m,k=1)(tvec)