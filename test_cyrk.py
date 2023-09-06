import numpy as np
from scipy.integrate import solve_ivp

c = np.asarray([(4 -  6**(1/2))/10, (4 +  6**(1/2))/10, 1])
a = np.asarray([[(88 - 7 * 6**(1/2))/360, (296 - 169 * 6**(1/2))/1800, (-2 + 3 * 6**(1/2))/225],
    [(296 + 169 * 6**(1/2))/1800, (88 + 7 * 6**(1/2))/360, (-2 - 3 * 6**(1/2))/225],
    [(16 -  6**(1/2))/36, (16 + 6**(1/2))/36, 1/9]])
b = np.asarray([(16 -  6**(1/2))/36, (16 +  6**(1/2))/36, 1/9])


def radau(ode, t0, y0, tf, h):
    t = t0
    y = y0
    while t < tf:
        k = []
        for i in range(len(c)):
            ti = t + c[i] * h
            yi = y + sum(a[i][j] * k[j] for j in range(i)) # must be solved implicitly
            k.append(ode(ti, yi))
        y = y + h * sum(b[i] * k[i] for i in range(len(b)))
        t = t + h
    return y

t0 = 0
tf = 1
h = 0.1
y0 = 1

def ode(t,y):
    return 1/np.sqrt(y)

solution = radau(ode, t0, y0, tf, h)

print(solve_ivp(ode,[t0,tf],[y0],"Radau",[tf])["y"])
print(solution)