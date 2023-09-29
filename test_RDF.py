import numpy as np
nP=200
zP=np.linspace(0,1,nP)
A=np.zeros((nP,nP))
C=np.zeros((nP,nP))
D=np.zeros((nP,nP))

# vec=np.asarray(range(0,nP*2,2))

h=1E-22
for i in range(nP):
    def RBF(zPj):
        return np.exp(-(epsilon*(zPj-zP[i]))**2)
         
    epsilon=0.4*nP
    for j in range(nP):
        # A[i,j]=np.exp(-(epsilon*(zP[j]-zP[i]))**2)
        A[i,j]=RBF(zP[j])
        # D[i,j]=-2*epsilon**2*(zP[j]-zP[i])*np.exp(-(epsilon*(zP[j]-zP[i]))**2)
        C[i,j]=np.imag(RBF(zP[j]+h*1j))/h
Ainv=np.linalg.inv(A)
S=C@Ainv
T=D@Ainv
print(S-T)
import matplotlib.pyplot as plt
[plt.plot(zP,A[i,:]) for i in range(0,nP,40)]
plt.show()