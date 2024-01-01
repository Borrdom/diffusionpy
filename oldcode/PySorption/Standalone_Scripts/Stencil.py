import numpy as np
def pointstencil(a):
    lower=a[0]-a[1]
    middle=-11/6*a[3:]+3*a[2:-1]-3/2*a[1:-2]+1/3*a[:-3]
    upper=a[-2]-a[-1]
    return b

a1=np.linspace(0,1,100)
a=a1**2
t=np.linspace(0,1,100)
lower=a[0]-a[1]
middle=-11/6*a[3:]+3*a[2:-1]-3/2*a[1:-2]+1/3*a[:-3]#4th order accuracy
upper=a[-2]-a[-1]
stencil=np.hstack((lower,middle,upper))