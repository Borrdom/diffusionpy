import numpy as np
import matplotlib.pyplot as plt
num=11
x1vec=np.linspace(0,1,num)
xspace=np.array([])
for i,val in enumerate(x1vec):
    x2=np.linspace(0,1-val,num-i)
    x1=val*np.ones_like(x2)
    x3=1-x2-x1
    x=np.asarray([x1,x2,x3])
    xspace=np.hstack((xspace,x)) if xspace.size else x
plt.plot(xspace[0,:],xspace[1,:],'kx')


dim=3
x1vec2=np.linspace(0,1,dim-2)

num=11
x1vec,step=np.linspace(0,1,num,retstep=True)

