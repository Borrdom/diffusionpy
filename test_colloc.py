import numpy as np

def broyden(fun,x):
    """solve this non-linear equation system for the nonbonded association sites XAi""" 
    tol=1E-8
    iter=50
    n=len(x)
    f=fun(x)
    h = tol
    J=np.ones((n,n)).astype(x.dtype)
    for i in range(n):
        dx = x*0
        dx[i] = h
        J[:,i] = (fun(x + dx)-f)/h
    for i in range(iter):
        if np.linalg.norm(f,2)<tol:
            return x 
        s=np.linalg.solve(J+np.eye(n)*tol,-1.*f)
        ff=fun(x+s)
        lamb=1.
        for j in range(iter):
            if np.linalg.norm(ff,2)*(1.-lamb*1E-1)<np.linalg.norm(f,2): break 
            lamb*=1./4.
            s*=lamb 
            ff=fun(x+s) 
        x+=s
        df=fun(x)-f
        J+=np.outer((df-np.dot(J,s)),s)/np.dot(s,s)
        f+=df
    return x


def orthogonal_collocation(n=4):
    z=(np.polynomial.legendre.leggauss(n+1)[0]+1)/2
    # z=np.asarray([0,0.1127,0.5,0.8873,1])
    A=np.zeros((n+1,n+1))
    C=np.zeros((n+1,n+1))
    D=np.zeros((n+1,n+1))
    for i in range(n+1):
        for j in range(n+1):
            A[i,j]=z[i]**j
            C[i,j]=j*z[i]**(j-1) if j>0 else 0
            D[i,j]=j*(j-1)*z[i]**(j-2) if j>1 else 0
    Ainv=np.linalg.inv(A)
    S=C@Ainv
    T=D@Ainv
    return z,T,S
def orthogonal_collocation_finite_elements(nP=4,nE=4):
    zE=np.linspace(0,1,nE+1)
    NS=np.polynomial.legendre.leggauss(nP)[0]
    NS=(NS-NS[0])/(NS[-1]-NS[0])
    n=(nP-1)*nE+1
    z=np.asarray([0])
    zvec=[]
    for i in range(nE):
        z0=zE[i]
        z8=zE[i+1]
        zP=NS*(z8-z0)+z0
        zvec.append(zP)
        z=np.hstack((z,zP[1:]))

    # z=np.asarray([0,0.1127,0.5,0.8873,1])
    A=np.zeros((nP,nP))
    C=np.zeros((nP,nP))
    D=np.zeros((nP,nP))
    for i in range(nP):
        for j in range(nP):
            A[i,j]=zP[i]**j
            C[i,j]=j*zP[i]**(j-1) if j>0 else 0
            D[i,j]=j*(j-1)*zP[i]**(j-2) if j>1 else 0
    Ainv=np.linalg.inv(A)
    S=C@Ainv
    T=D@Ainv

    return z,T,S



nP=4
nE=4
n=(nP-1)*nE+1
# z,T,S=orthogonal_collocation(n)
z,T,S=orthogonal_collocation_finite_elements(nP,nE)
def fun_(c,dcdz,d2cd2z):
    Pe=1
    Da=1
    f=1/Pe*d2cd2z-dcdz-Da*c**2
    return f 
def fun(c):
    f=np.zeros_like(c)
    Pe=1
    Da=1

    for i in range(nE):
        # f[i*nP] # übergang?
        P0=i*(nP-1)
        P8=(i+1)*(nP-1)+1
        cP=c[P0:P8]
        
        d2cd2z=T@cP
        dcdz=S@cP

        f[P0:P8]=fun_(cP,dcdz,d2cd2z)
        if i>0: f[(nP-1)*i]=dcdz_[-1]-dcdz[0] 
        if i==0: f[0]=dcdz[0]-Pe*(c[0]-1) # RB1 
        if i==nE: f[-1]=dcdz[-1]# RB2
        dcdz_=dcdz
    return f
# def fun(c):
#     Pe=1
#     Da=1
#     f=np.zeros_like(c)
#     for i in range(n+1):
#         f[i]=1/Pe*np.dot(T[i,:],c)-np.dot(S[i,:],c)-Da*c[i]**2
#     f[0]=np.dot(S[0,:],c)-Pe*(c[0]-1)
#     f[-1]=np.dot(S[-1,:],c)
#     return f
c0=np.ones(n)
c=broyden(fun,c0)
import matplotlib.pyplot as plt
plt.plot(z,c)
plt.show()
print(c)