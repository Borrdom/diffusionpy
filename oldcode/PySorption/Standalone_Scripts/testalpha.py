import numpy as np
Xa=np.linspace(0,1,100)
Xa0=0.7
alpha=(Xa0-Xa)/(1-Xa)
import matplotlib.pyplot as plt
#plt.plot(Xa,alpha)
#plt.show()
XWeq=0.3
Xaeq=0.05
alphaeq=(Xa-Xaeq)/(1-Xaeq)
XWL=XWeq/(1-Xaeq)**2*(1-Xa)**2
XW=(1-alphaeq)*XWeq
XW2=(1-Xa)/(1-Xaeq)*XWeq
ww2=XW2/(1+XW2)
ww=XW/(1+XW)
wwL=XWL/(1+XWL)
#plt.plot(Xa,XWL,"k-")
#plt.plot(Xa,XW,'b-')
#plt.plot(Xa,XW2,'b--')
#plt.show()
plt.plot(Xa,wwL,"k-")
plt.plot(Xa,ww,'b-')
plt.plot(Xa,ww2,'b--')
plt.show()
