
import numpy as np
import matplotlib.pyplot as plt
k=1
n=np.linspace(0,1,100)
t=np.linspace(0,1,10000)
m=lambda t,n:k*t**n
fig,ax=plt.subplots(2,1)
ax[0].set_xlabel("t^0.5 [s]")
ax[0].set_ylabel("m [-]")
ax[1].set_xlabel("t [s]")
ax[1].set_ylabel("m [-]")
[ax[0].plot(np.sqrt(t),m(t,val)) for i,val in enumerate(n)]
[ax[1].plot(t,m(t,val)) for i,val in enumerate(n)]