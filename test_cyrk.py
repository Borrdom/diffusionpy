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


import numpy as np
from matplotlib import pyplot as plt

x0 = 0
y0 = 2
xf = 1
n = 11
deltax = (xf-x0)/(n-1)
x = np.linspace(x0,xf,n)
def f(x,y):
	return y-x

y = np.zeros([n])
y[0] = y0
py = np.zeros([n])
for i in range(0,4):
	py[i] = None

for i in range(1,4):
	k1 = deltax*f(x[i-1],y0)
	k2 = deltax*f(x[i-1]+deltax/2,y0+k1/2)
	k3 = deltax*f(x[i-1]+deltax/2,y0+k2/2)
	k4 = deltax*f(x[i-1]+deltax,y0+k3)
	y[i] =  y0 + (k1 + 2*k2 + 2*k3 + k4)/6
	y0 = y[i]

for i in range(4,n):
	py[i] = deltax/24*(55*f(x[i-1],y[i-1]) - 59*f(x[i-2],y[i-2]) + 37*f(x[i-3],y[i-3]) - 9*f(x[i-4],y[i-4]) )  + y[i-1] 
	y[i] = deltax/24*( 9*f(x[i],py[i]) + 19*f(x[i-1],y[i-1]) - 5*f(x[i-2],y[i-2]) + f(x[i-3],y[i-3]) ) + y[i-1]

print("x_n\t   py_n\t           y_n")
for i in range(n):
	print (x[i],"\t",format(py[i],'6f'),"\t",format(y[i],'6f'))

plt.plot(x,y,'o')
plt.xlabel("Value of x")
plt.ylabel("Value of y")
plt.title("Approximation Solution with Adams-Bashforth-Moulton Method")
plt.show()