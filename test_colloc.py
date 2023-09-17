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
    NS=np.polynomial.chebyshev.chebpts1(nP)
    # NS=np.polynomial.legendre.leggauss(nP)[0]
    NS=(NS-NS[0])/(NS[-1]-NS[0])
    print(NS.tolist())
    # NS=np.asarray([0,0.1127,0.5,0.8873,1])
    NS=np.linspace(0,1,nP) # equidistant collocation
    n=(nP-1)*nE+1
    z=np.asarray([0])
    zvec=[]
    for i in range(nE):
        z0=zE[i]
        z8=zE[i+1]
        zP=NS*(z8-z0)+z0
        zvec.append(zP)
        z=np.hstack((z,zP[1:]))

    
    A=np.zeros((nP,nP))
    C=np.zeros((nP,nP))
    D=np.zeros((nP,nP))
    # vec=np.asarray(range(0,nP*2,2))
    for i in range(nP):
        for j in range(nP):
            # A[i,j]=NS[i]**j
            # C[i,j]=j*NS[i]**(j-1) if j>0 else 0
            # D[i,j]=j*(j-1)*NS[i]**(j-2) if j>1 else 0
    # for i,vali in enumerate(vec):
    #     for j,valj in enumerate(vec):
    #         A[i,j]=zP[i]**valj
    #         C[i,j]=valj*zP[i]**(valj-1) if j>0 else 0
    #         D[i,j]=valj*(valj-1)*zP[i]**(valj-2) #if j>1 else 0
            A[i,j]=np.exp(zP[i]*j)
            C[i,j]=j*np.exp(zP[i]*j)
            D[i,j]=j**2*np.exp(-zP[i]*j)
    Ainv=np.linalg.inv(A)
    S=C@Ainv
    T=D@Ainv

    return z,T,S

import matplotlib.pyplot as plt
nshow=4
fig,ax=plt.subplots(1,nshow)


for w in range(nshow):
    nP=2
    nE=1+w
    n=(nP-1)*nE+1
    # z,T,S=orthogonal_collocation(n)
    z,T,S=orthogonal_collocation_finite_elements(nP,nE) 
    print(S.tolist())
    def fun_(c,dcdz,d2cd2z,Pe,Da):
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
            
            dcdz=(S@cP)
            d2cd2z=(S@dcdz)
            # d2cd2z=T@cP
            if i>0: dcdz[0]=dcdz_[-1]
            f[P0:P8]=fun_(cP,dcdz,d2cd2z,Pe,Da)
            
            # if i>0: f[(nP-1)*i]=dcdz_[-1]-dcdz[0] 
            if i==0: f[0]=dcdz[0]-Pe*(c[0]-1) # RB1 
            if i==(nE-1): f[-1]=dcdz[-1]# RB2
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
    
    
    ax[w].plot(z,c,"o-")
    zE=np.linspace(0,1,nE+1)
    for i in range(nE):
        ax[w].plot([zE[i],zE[i]],[0.59,0.7311],"k")
        ax[w].plot([zE[i+1],zE[i+1]],[0.59,0.7311],"k")
plt.show()
    # print(c)
S=np.asarray([[-9.183300132670368, 13.95870750793668, -7.32835571279449, 3.5529483375281643, -0.9999999999999911],  # Collocation matrix calculated for 5 collocation points
[-1.7403293111129658, -1.3744088030005464, 4.3546484316145255, -1.682881139058434, 0.44297082155741346], 
[0.5458250331675965, -2.601439792602136, 3.552713678800501e-15, 2.6014397926021253, -0.5458250331675947], 
[-0.4429708215574084, 1.6828811390584317, -4.354648431614518, 1.3744088030005417, 1.7403293111129672], 
[0.9999999999999858, -3.552948337528136, 7.328355712794462, -13.95870750793668, 9.183300132670368]])    

S1=np.asarray([[-9.183300132670368, 13.95870750793668, -7.32835571279449, 3.5529483375281643, -0.9999999999999911], [-1.7403293111129658, -1.3744088030005464, 4.3546484316145255, -1.682881139058434, 0.44297082155741346], [0.5458250331675965, -2.601439792602136, 3.552713678800501e-15, 2.6014397926021253, -0.5458250331675947], [-0.4429708215574084, 1.6828811390584317, -4.354648431614518, 1.3744088030005417, 1.7403293111129672], [0.9999999999999858, -3.552948337528136, 7.328355712794462, -13.95870750793668, 9.183300132670368]])
S2=np.asarray([[-18.366600265340224, 27.91741501587171, -14.656711425586707, 7.105896675054453, -1.9999999999994031], [-3.480658622226039, -2.748817606001179, 8.709296863229428, -3.3657622781170686, 0.885941643114839], [1.0916500663349318, -5.202879585203583, -7.815970093361102e-14, 5.202879585204428, -1.0916500663353403], [-0.8859416431147413, 3.3657622781166485, -8.709296863228746, 2.748817606000955, 3.4806586222259464], [2.0000000000002274, -7.105896675056556, 14.656711425589947, -27.917415015873303, 18.366600265340708]])
S3=np.asarray([[-27.549900397995966, 41.876122523767556, -21.98506713833475, 10.658845012549834, -2.9999999999884817], [-5.220987933338487, -4.123226409005926, 13.063945294846087, -5.048643417178575, 1.328912464673726], [1.6374750994997316, -7.804319377800917, -1.1498676291743795e-11, 7.804319377812084, -1.6374750995057181], [-1.3289124646725352, 5.04864341717289, -13.063945294847636, 4.123226408999917, 5.2209879333387], [2.9999999999936335, -10.658845012572783, 21.985067138366503, -41.87612252380313, 27.54990039800896]])
S4=np.asarray([[-36.73320053070856, 55.834830031822094, -29.313422851292216, 14.21179335018951, -4.000000000030553], [-6.96131724444397, -5.4976352120264345, 17.4185937264693, -6.731524556254808, 1.7718832862322151], [2.1833001326761554, -10.405759170429405, 7.673861546209082e-12, 10.405759170396891, -2.183300132668304], [-1.7718832862290812, 6.731524556227857, -17.418593726464845, 5.4976352120031216, 6.9613172444471045], [4.000000000016371, -14.211793350172229, 29.313422851228097, -55.83483003179572, 36.733200530696195]])
# S=np.asarray([
    # [-25./3., 16., -12., 16./3., -1.],
    # [-1., -10./3., 6., -2., 1./3.],
    # [1/3, -8./3., 0., 8./3, -1./3.], 
    # [-1./3., 2., -6., 10./3., 1.],
    # [1., -16./3., 12., -16., 25./3.]]) #equidistant collocation so that  the equations are inependent of the length between the points