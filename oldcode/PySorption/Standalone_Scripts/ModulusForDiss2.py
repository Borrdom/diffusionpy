#import numpy as np
from scipy.integrate import odeint
import autograd.numpy as np
from autograd import grad 
import matplotlib.pyplot as plt
from PySorption.DiffusionsmodellPlots import Plotlist



nArm=2
logt=np.linspace(-2,3,300)
t=10**logt
epsilon0=1
taudiff=10**2
epsilon_fun1=lambda t : (1-np.exp(-t/taudiff))*epsilon0
#epsilon_fun1=lambda t : np.fmin(epsilon0*t,1)
#epsilon_fun1=lambda t : np.fmin(epsilon0*t,1)
epsilon_dot_fun1= grad(epsilon_fun1)
#epsilon_dot_fun1=lambda t: 0
#epsilon_dot_fun2= grad(epsilon_fun2)
tauJ=np.asarray([10**-2,10**2])
EJ=np.asarray([1E9,1E5])

etaJ=tauJ*EJ
#pK=tauJ
#qK=EJ

#Num=sum([0+val*(omega*1j)**k for k,val in enumerate(pK)])
#Den=sum([0+val*(omega*1j)**k for k,val in enumerate(qK)])

#E=Num/Den
#E=sum([Eeq+(omega*1j)/(1+val*omega*1j)*etaJ[k] for k,val in enumerate(tauJ)])
#sigma1=tauJ*epsilon0*EJ*(1-np.exp(-t/tauJ))
#sigma2=epsilon0*EJ*(np.exp(-t/tauJ))
def MaxwellElement(y,t,EJ,tauJ,epsilon_dot_fun,nArm):
    epsilon=np.asarray([y[i] for i in range(nArm)])
    sigma=np.asarray([y[nArm+i] for i in range(nArm)])
    epsilon_dot=np.asarray([epsilon_dot_fun(t) for i in range(nArm)])
    sigma_dot=EJ*epsilon_dot-sigma/tauJ
    y_dot=np.hstack((epsilon_dot,sigma_dot))
    return y_dot
    
y0=np.hstack((np.zeros(nArm),np.ones(nArm)))
sol1=odeint(MaxwellElement,y0,t,args=(EJ,tauJ,epsilon_dot_fun1,nArm))

epsilon1=sol1[:,0]
epsilon2=sol1[:,1]
sigma1=sol1[:,2]
sigma2=sol1[:,3]

E1=(sigma1+sigma2)/epsilon1


#storage=E.real
#loss=E.imag
#tandelta=storage/loss
fig,ax=plt.subplots()
ax.plot(t,sigma2,'r-')

fig1,ax1=plt.subplots()
ax1.plot(t,epsilon2,'r-')
#ax.plot(logt,sigma2,'b-')
#ax.plot(logomega,loss,'b-')
#ax.plot(logomega,tandelta,'g-')
x=[t,t]
y=[epsilon1,E1]
xs=x
ys=y
#Plotlist(x,y,xs,ys,xe=0,ye=0,Origin=True,filename="Emoddis",xla="t/s",yla="E/Pa")
plt.show()