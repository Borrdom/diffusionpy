import numpy as np

ww=np.linspace(0,1)
D0=1E-13
vH2O=0.1
vPol=0.5
v=vH2O+vPol
a=10
b=1
zeta=-a*ww+b

D=D0*np.exp(-(vH2O*ww+zeta*vPol*(1-ww))/v)
import matplotlib.pyplot as plt

plt.plot(ww,D)
plt.show()
