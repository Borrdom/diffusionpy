import casadi as cs
import numpy as np
import matplotlib.pyplot as plt
class Complex(object):
    def __init__(self, re,im):
        self.re = re
        self.im = im

    # adding two objects
    def __add__(self, o):
        return Complex(self.re + o.re,+self.im + o.im)
    def __sub__(self, o):
        return Complex(self.re - o.re,+self.im - o.im)
    def __mul__(self, o):
        return Complex((self.re * o.re)-(self.im*o.im),(self.re*o.im)+(self.im * o.re))
    def __truediv__(self, o):
        den=o.re * o.re + o.im * o.im
        num1=self.re * o.re+self.im * o.im
        num2=self.im * o.re-self.re * o.im
        return Complex(num1/den,+num2/den)
    def __pow__(self,o):
        abs=(self.re*self.re+self.im*self.im)**0.5
        abso=abs**o
        theta=cs.arctan(self.im/self.re)

        return Complex(abso*cs.cos(theta*o),abso*cs.sin(theta*o))

 
if __name__=="__main__":
    A=Complex(5,1)
    B=Complex(1,2)
    a=5+1j
    b=1+2j


    C=A*B
    D=A/B
    c=a*b
    d=a/b
    E=A**3
    e=a**3
    print(C.re-c.real)
    print(C.im-c.imag)

    print(D.re-d.real)
    print(D.im-d.imag)

    print(E.re-e.real)
    print(E.im-e.imag)
    nz=1000
    X=cs.MX.sym("X",nz)
    L=1
    dz=L/nz
    def CDF(a):
        diff=(a[1:]-a[:-1])/dz
        #diff=Complex(a[1:],0).im
        return diff 
        
    #dXdt=cs.vertcat(CDF(cs.vertcat(0,CDF(X))),0)
    m=0
    for i in range(nz-1):
        Xco=Complex(X[i+1],dz)
        m=cs.vertcat(m,Xco.im/dz)

    dXdt=[]
    for i in range(nz-1):
        mco=Complex(m[i+1],dz)
        dXdt=cs.vertcat(dXdt,mco.im/dz)
    dXdt=cs.vertcat(dXdt,0)

    ode_def = {'x':X,'ode':dXdt}
    nt=100
    ts=np.linspace(0,10,nt)
    opts = {}
    opts ['abstol'] = 1e-6
    opts['reltol'] = 1e-6
    opts['grid'] = ts
    opts['output_t0'] = True

    F = cs.integrator('sim', 'idas', ode_def, opts)
    x0=np.zeros(nz)
    x0[-1]=1
    sol = F(x0=x0)
    solx = np.reshape(sol['xf'].full(),(nt,nz))
    [plt.plot(ts,solx[:,i]) for i,val in enumerate(solx[0,:])]
    print(solx)
    plt.show()