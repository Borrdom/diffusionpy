import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

fig = plt.figure()
ax = fig.add_subplot(111)


# Plot the surface.



u=np.linspace(200,400,100)
t=np.linspace(50,500,100)
U,T=np.meshgrid(u,t)

factor=(1-0.12*np.exp(-3*U/T))**3
factor1=(1-0.12*np.exp(-3*250/298.15))**3
factor2=(1-0.12*np.exp(-3*350/298.15))**3
turf=ax.contourf(U, T, factor)
ax.plot(250,298.15,'rx')
ax.plot(350,298.15,'gx')
ax.set_xlabel('uij/kb/-')
ax.set_ylabel('T/K')
plt.colorbar(turf)

plt.show()
