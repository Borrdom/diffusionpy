
import numpy as np
from scipy.integrate import odeint

def model(y,t):
    k=0.0003
    dydt=-k*y*t
    return dydt
#model= lambda y,t:-k*t*y
y0=0.0001
t0=0
tend=200
nt=50
t=np.linspace(t0,tend,nt)
y=odeint(model,y0,t)
