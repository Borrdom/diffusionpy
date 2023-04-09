import numpy as np
from numba import njit, config
import matplotlib.pyplot as plt
h=1E-26
kB = 1.380649e-23
NA = 6.0221407e23
R=kB*NA


# config.DISABLE_JIT = True

#@njit(optional(complex128)(float64,optional(complex128),optional(complex128[:,:]),float64[:,:],float64[:,:],float64[:,:],float64[:,:],float64[:,:],float64[:,:],float64[:,:]),fastmath=True)

#@njit(cache=True)
@njit(['float64(float64,float64,float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1])',
        'complex128(float64,complex128,float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1])',
        'complex128(float64,float64,complex128[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1])'],
cache=True)
def ares(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N):

    def npaddouter(a): return a.reshape(len(a),1)+a

    def broyden(fun,x,p1,p2,p3,tol=1E-8,iter=50):
        n=len(x)
        f=fun(x,p1,p2,p3)#.reshape(n)
        h = tol
        J=np.zeros((n,n))

        for i in range(n):
            dx = x*0
            dx[i] = h #* 1j     
            #J[i,:] = (fun(x + dx,p1,p2,p3).reshape(n)-f)/h
            J[i,:] = (fun(x + dx,p1,p2,p3)-f)/h
        lamb=f*0+1.
        for i in range(iter):
            if np.linalg.norm(f,n)<tol:
                return x 
            s=np.linalg.solve(J+np.eye(n)*1E-6,-1.*f)
            ff=fun(x+s,p1,p2,p3)#.reshape(n)
            for j in range(n):
                if np.abs(ff[j])>np.abs(f[j]):
                    lamb[j]*=1./4.
                    s[j]*=lamb[j]    
            x+=s#.reshape(n,1)
            #df=fun(x,p1,p2,p3).reshape(n)-f
            df=fun(x,p1,p2,p3)-f
            J+=np.outer((df-np.dot(J,s)),s)/np.dot(s,s)
            f+=df
            #plt.plot(i,np.log10(np.linalg.norm(f,n)),'kx')
        return x

    npolyI=7
    npolycoef=3

    a0=np.ones((npolycoef,npolyI))
    b0=np.ones((npolycoef,npolyI))
    a0[:,0]=np.asarray([0.9105631445,-0.3084016918,-0.0906148351])
    a0[:,1]=np.asarray([0.6361281449,0.1860531159,0.4527842806])
    a0[:,2]=np.asarray([2.6861347891,-2.5030047259,0.5962700728])
    a0[:,3]=np.asarray([-26.547362491,21.419793629,-1.7241829131])
    a0[:,4]=np.asarray([97.759208784,-65.255885330,-4.1302112531])
    a0[:,5]=np.asarray([-159.59154087,83.318680481,13.776631870])
    a0[:,6]=np.asarray([91.297774084,-33.746922930,-8.6728470368])

    b0[:,0]=np.asarray([0.7240946941,-0.5755498075,0.0976883116])
    b0[:,1]=np.asarray([2.2382791861,0.6995095521,-0.2557574982])
    b0[:,2]=np.asarray([-4.0025849485,3.8925673390,-9.1558561530])
    b0[:,3]=np.asarray([-21.003576815,-17.215471648,20.642075974])
    b0[:,4]=np.asarray([26.855641363,192.67226447,-38.804430052])
    b0[:,5]=np.asarray([206.55133841,-161.82646165,93.626774077])
    b0[:,6]=np.asarray([-355.60235612,-165.20769346,-29.666905585])
    ncomp=len(mi)
    ntype=2

    #Initializekij
    kij=np.zeros((ncomp,ncomp))
    lij=np.zeros((ncomp,ncomp))

    # Hard Chain Contribution
    mibar=np.sum(xi*mi)
    di=sigi*(1.-0.12*np.exp(-3*ui/T))
    rho=6/np.pi*eta*(np.sum(mi*xi.real*di**3))**-1 #mi*xi=phi*mbar und dann mbar*rho=rhophi
    
    z3=np.pi/6*rho*np.sum(mi*xi*di**3)
    z2=np.pi/6*rho*np.sum(mi*xi*di**2)
    z1=np.pi/6*rho*np.sum(mi*xi*di**1)
    z0=np.pi/6*rho*np.sum(mi*xi*di**0)
    #Mixing rules

    dij=npaddouter(di)
    sigij=npaddouter(sigi)/2
    #uij=(ui@ui.T)**0.5*(1.-kij)
    uij=np.outer(ui,ui)**0.5*(1.-kij)
    #kapij=(kapi@kapi.T)**0.5*((sigi@sigi.T)**0.5/sigij)**3
    kapij=np.outer(kapi,kapi)**0.5*(np.outer(sigi,sigi)**0.5/sigij)**3
    epsAiBj=npaddouter(epsAiBi)/2.
    #gij=1/(1-z3)+(di@di.T/dij)*3*z2/(1-z3)**2+((di@di.T)/(dij))**2*2*z2**2/(1-z3)**3
    gij=1/(1-z3)+(np.outer(di,di)/dij)*3*z2/(1-z3)**2+(np.outer(di,di)/(dij))**2*2*z2**2/(1-z3)**3
    fhs=1/z0*(3*z1*z2/(1-z3)+z2**3/(z3*(1-z3)**2)+(z2**3/z3**2-z0)*np.log(1-z3))
    #diaggij=np.diag(gij).reshape(ncomp,1)
    diaggij=np.diag(gij)#.reshape(ncomp,1)

    fhc=mibar*fhs-np.sum(xi*(mi-1)*np.log(diaggij))

    # Dispersion Contribution
    a=a0[0,:]+(mibar-1)/mibar*a0[1,:]+(mibar-1)/mibar*(mibar-2)/mibar*a0[2,:]
    b=b0[0,:]+(mibar-1)/mibar*b0[1,:]+(mibar-1)/mibar*(mibar-2)/mibar*b0[2,:]

    I1=a[0]+a[1]*z3+a[2]*z3**2+a[3]*z3**3+a[4]*z3**4+a[5]*z3**5+a[6]*z3**6
    I2=b[0]+b[1]*z3+b[2]*z3**2+b[3]*z3**3+b[4]*z3**4+b[5]*z3**5+b[6]*z3**6


    C1=(1+mibar*(8*z3-2*z3**2)/(1-z3)**4+(1-mibar)*(20*z3-27*(z3**2)+12*(z3**3)-2*z3**4)/((1-z3)*(2-z3))**2)**-1

    #m2es3mat=(xi@xi.T)*(mi@mi.T)*(uij/T)*sigij**3
    m2es3mat=np.outer(xi,xi)*np.outer(mi,mi)*(uij/T)*sigij**3

    m2e2s3mat=m2es3mat*(uij/T)
    m2es3,m2e2s3=np.sum(m2es3mat),np.sum(m2e2s3mat)

    fdisp=-2*np.pi*rho*I1*m2es3-np.pi*rho*mibar*C1*I2*m2e2s3
    

    # Association Contribution
    deltAiBj=gij*kapij*sigij**3*(np.exp(epsAiBj/T)-1.)
    rhoi=rho.real*xi.real
    #def XAi_eq(XAi,rhoi,N,deltAiBj): return XAi-((1+np.sum(rhoi*XAi*N*deltAiBj.T,axis=0))**-1).reshape(ncomp,1)
    def XAi_eq(XAi,rhoi,N,deltAiBj): return XAi-((1+np.sum(rhoi*XAi*N*deltAiBj.T,axis=1))**-1)
    #XAiself=((-1+np.sqrt(1+4*rho.real*np.diag(deltAiBj.real)))/(2*rho.real*np.diag(deltAiBj.real))).reshape(ncomp,1)
    XAiself=((-1.+np.sqrt(1.+4.*rho.real*np.diag(deltAiBj.real)))/(2.*rho.real*np.diag(deltAiBj.real)))
    XAiself[np.isnan(XAiself)]=1.
    
    #XAiself=((-1+np.sqrt(1+8*rhoi.flatten()*np.diag(deltAiBj.real)))/(4*rhoi.flatten()*np.diag(deltAiBj.real))).reshape(ncomp,1)
    #XAiself=np.ones_like(xi.real)*0.5
    XAi=broyden(XAi_eq,XAiself,rhoi,N,deltAiBj.real)#.reshape(ncomp,1)

    tm4=np.sum((np.log(XAi)-XAi+1)*N*ntype*xi)
    tmp4=np.sum(np.outer(XAi,XAi)*deltAiBj*np.outer(N,N)*np.outer(xi,xi)*ntype)
    q=tm4-rho/2*tmp4
    
    fassoc=np.sum((np.log(XAi)-1/2*XAi+1/2)*N*ntype*xi)
    

    
    return fhc+fdisp+q

@njit('float64(float64,float64,float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1])',cache=True)
def Z(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N):
    return 1.+eta*(ares(T,eta+h*1.j,xi,mi,sigi,ui,epsAiBi,kapi,N).imag)/h

@njit('float64[::1](float64,float64,float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1])',cache=True)
def lnphii(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N):
    ncomp=len(xi)
    #dadx=np.ones((ncomp,1))
    dadx=np.ones(ncomp)
    a1=ares(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N)
    Z1=Z(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N)
    hvec=np.eye(ncomp)*h*1.j
    for i in range(ncomp):
        xi_=xi+hvec[i,:]#.reshape(ncomp,1)
        dadx[i]=ares(T,eta,xi_,mi,sigi,ui,epsAiBi,kapi,N).imag/h
    lnphi = a1 + (Z1-1) + dadx - np.dot(xi, dadx)#  -np.log(Z1)
    return lnphi

# #@njit(cache=True)

# def XAi_obj(T,eta,xi,XAi,mi,sigi,ui,epsAiBi,kapi,N):
#     di=sigi*(1-0.12*np.exp(-3*ui/T))
#     rho=6/np.pi*eta*(np.sum(mi*xi*di**3))**-1 #mi*xi=phi*mbar und dann mbar*rho=rhophi    
#     z3=np.pi/6*rho*np.sum(mi*xi*di**3)
#     z2=np.pi/6*rho*np.sum(mi*xi*di**2)
#     ncomp=len(xi)
#     #Mixing rules
#     dij=matmultadd(di)
#     sigij=matmultadd(sigi)/2
#     kapij=(kapi@kapi.T)**0.5*((sigi@sigi.T)**0.5/sigij)**3
#     epsAiBj=matmultadd(epsAiBi)/2
#     gij=1/(1-z3)+(di@di.T/dij)*3*z2/(1-z3)**2+((di@di.T)/(dij))**2*2*z2**2/(1-z3)**3    
#     # Association Contribution
#     deltAiBj=gij*kapij*sigij**3*(np.exp(epsAiBj/T)-1.)
#     XAi_eq=(1+np.sum(rho*xi*XAi*N*deltAiBj.T,axis=0))**-1
#     return XAi_eq.reshape(ncomp,1)-XAi
#@njit(cache=True)
def Z_obj(p,T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N):
    di=sigi*(1.-0.12*np.exp(-3*ui/T))
    rho=6/np.pi*eta*(np.sum(mi*xi*di**3))**-1
    rhobar=rho*(10.**10)**3
    Zp=p/(rhobar*kB*T)
    return (Zp-Z(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N))

# def lnphi_p(p,T,xi,mi,sigi,ui,epsAiBi,kapi,N):
#     def TP_obj(eta): return Z_obj(p,T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N)
#     eta0=0.4
#     xopt=root(TP_obj,eta0, method='lm')["x"]
#     eta_n=xopt[0]
#     return lnphii(T,eta_n,xi,mi,sigi,ui,epsAiBi,kapi,N)

@njit('float64[::1](float64,float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1],float64[::1])',cache=True)
def SAFTSAC(T,vpure,xi,mi,sigi,ui,epsAiBi,kapi,N):
    vmol=np.sum(vpure*xi)
    vpfrac=vpure/vmol
    #rho=6/np.pi*eta*(np.sum(mi*xi*di**3))**-1
    di=sigi*(1.-0.12*np.exp(-3*ui/T))
    eta=np.pi/6*np.sum(mi*xi*di**3)/vmol/(10.**10)**3*NA
    etapure=np.pi/6*mi*di**3/vpure/(10.**10)**3*NA
    lngammaid=np.log(vpfrac)+1-vpfrac
    arespures=np.asarray([ares(T,val,np.asarray([1.]),np.asarray([mi[i]]),np.asarray([sigi[i]]),np.asarray([ui[i]]),np.asarray([epsAiBi[i]]),np.asarray([kapi[i]]),np.asarray([N[i]])) for i,val in enumerate(etapure)])
    lngammares=lnphii(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N)-arespures
    p=Z(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N)/vmol*(R*T)
    lngammap=vpure/(R*T)*(p-(R*T/vmol))
    return lngammaid+lngammares-lngammap




#Test call, so it compiles directly when it is imported
T=298.15
p=1E5
npoint=1
#Water
#XAiself =0.16 und 0.6
mi=np.asarray([1.20469,2.38269789])
sigi=np.asarray([2.797059952,3.1771])
ui=np.asarray([353.95,198.24])
epsAiBi=np.asarray([2425.67,2653.4])
kapi=np.asarray([0.04509,0.032384])
N=np.asarray([1.,1.])

#
# mi=np.asarray([3.057603049,2.53030029])
# sigi=np.asarray([3.7983,3.8499])
# ui=np.asarray([236.77,278.11])
# epsAiBi=np.asarray([0.,0.])
# kapi=np.asarray([0.,0.])
# N=np.asarray([0.,0.])

xi=np.asarray([0.5,0.5])
x1=np.linspace(0,1,11)
x2=1-x1
xi=np.vstack((x1,x2))
vpure=np.asarray([(996.9651701/18.015*1000)**-1,(779.7859022/46.069*1000.)**-1])
#vpure=np.asarray([(649.6531684/86.177*1000)**-1,(762.4382449/84.147*1000.)**-1])


#lngammai=SAFTSAC(T,vpure,xi,mi,sigi,ui,epsAiBi,kapi,N).flatten()
lngammai=np.asarray([SAFTSAC(T,vpure,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N).flatten() for i,val in enumerate(xi[0,:])])
print(lngammai)