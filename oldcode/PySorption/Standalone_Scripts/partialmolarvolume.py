import numpy as np
import matplotlib.pyplot as plt

ww=np.linspace(0,1,10000)
wp=1-ww
rho0p=1150

rho0w=1000
rho=rho0p/wp
v=1/rho
vid=ww*1/rho0w+wp*1/rho0p

#plt.plot(ww,rho,"o")


plt.plot(ww,v,"k-")
plt.plot(ww,vid,"k--")
plt.show()
