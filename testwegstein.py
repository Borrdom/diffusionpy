import numpy as np
from numba import njit
import matplotlib.pyplot as plt

@njit
def wegstein(fun,x, p1):
    """Solving via wegsteins method"""
    tol=1E-8
    maxiter=50
    f = fun(x, p1)
    xx=f    
    dx = xx - x
    ff = fun(xx, p1)
    df=ff-f
    for i in range(maxiter):
        e=np.linalg.norm(df,2)
        print(f"iter {i}: ||F|| = {e}")
        if e<tol: 
            return xx 
        a = df/dx
        q= a/(a-1)
        x=xx
        xx = q * xx + (1-q) * ff
        f=ff
        ff = fun(xx, p1)    
        df=ff-f
        dx=xx-x
    return xx

@njit
def powell_fun(x,p1):
    f=np.zeros_like(x)
    f[0]=np.sinh(p1*x[0])
    f[1]=np.sinh(p1*x[1])
    return f

x0=np.asarray([1.,1.])
p1=1.2
f1=powell_fun(x0,p1)

print(f1)


fsol=wegstein(powell_fun,x0,p1)

print(fsol)
print(powell_fun(fsol,p1))
print(powell_fun(powell_fun(fsol,p1),p1))
plt.show()