import numpy as np
import matplotlib.pyplot as plt
nz=100
z=np.linspace(-0.5,0.5,nz)

rhov=1
rhol=997

b=10
def rho_(x):
    phi=(np.tanh(x*b)+1)/2
    return phi*rhol+(1-phi)*rhov

drho_dz=lambda x: rho_(x+1j*1E-22).imag/1E-22
rho=rho_(z)
plt.plot(z,rho)
plt.plot(z,drho_dz(z))

plt.show()



