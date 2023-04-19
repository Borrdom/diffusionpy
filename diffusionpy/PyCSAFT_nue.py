import numpy as np
from numba import njit
@njit(cache=True)
def ares(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB):
    def npaddouter(a): 
        return a.reshape(len(a),1)+a
    def wertheimiter(fun,x,p1,p2,p3,tol=1E-8,iter=50):
        n=len(x)
        f=fun(x,p1,p2,p3)
        h = tol
        J=np.ones((n,n)).astype(x.dtype)
        for i in range(n):
            dx = x*0
            dx[i] = h
            J[:,i] = (fun(x + dx,p1,p2,p3)-f)/h
        for i in range(iter):
            if np.linalg.norm(f,2)<tol:
                return x 
            s=np.linalg.solve(J+np.eye(n)*tol,-1.*f)
            ff=fun(x+s,p1,p2,p3)
            lamb=1.
            for j in range(iter):
                if np.linalg.norm(ff,2)*(1.-lamb*1E-1)<np.linalg.norm(f,2): break 
                lamb*=1./4.
                s*=lamb 
                ff=fun(x+s,p1,p2,p3) 
            x+=s
            df=fun(x,p1,p2,p3)-f
            J+=np.outer((df-np.dot(J,s)),s)/np.dot(s,s)
            f+=df
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
    kij=np.zeros((ncomp,ncomp)) if kij is None else kij
    kijAB=np.zeros((ncomp,ncomp)) if kijAB is None else kijAB

    #Mixing rules
    mibar=np.sum(xi*mi)
    di=sigi*(1.-0.12*np.exp(-3*ui/T))
    rho=6/np.pi*eta*(np.sum(mi*xi.real*di**3))**-1 
    dij=npaddouter(di)
    sigij=npaddouter(sigi)/2.
    uij=np.outer(ui,ui)**0.5*(1.-kij)
    kapij=np.outer(kapi,kapi)**0.5*(np.outer(sigi,sigi)**0.5/sigij)**3
    epsAiBj=npaddouter(epsAiBi)/2.*(1.-kijAB)

    # Hard-Chain Contribution
    z3=np.pi/6.*rho*np.sum(mi*xi*di**3)
    z2=np.pi/6.*rho*np.sum(mi*xi*di**2)
    z1=np.pi/6.*rho*np.sum(mi*xi*di**1)
    z0=np.pi/6.*rho*np.sum(mi*xi*di**0)   
    dijmat=np.outer(di,di)/dij
    diimat=np.diag(dijmat)
    gij=1./(1.-z3)+(dijmat)*3.*z2/(1.-z3)**2+(dijmat)**2*2*z2**2/(1.-z3)**3
    fhs=1./z0*(3.*z1*z2/(1-z3)+z2**3/(z3*(1.-z3)**2)+(z2**3/z3**2-z0)*np.log(1.-z3))
    gii=np.diag(gij)
    fhc=mibar*fhs-np.sum(xi*(mi-1)*np.log(gii))

    # Dispersion Contribution
    a=a0[0,:]+(mibar-1)/mibar*a0[1,:]+(mibar-1)/mibar*(mibar-2)/mibar*a0[2,:]
    b=b0[0,:]+(mibar-1)/mibar*b0[1,:]+(mibar-1)/mibar*(mibar-2)/mibar*b0[2,:]
    I1=a[0]+a[1]*z3+a[2]*z3**2+a[3]*z3**3+a[4]*z3**4+a[5]*z3**5+a[6]*z3**6
    I2=b[0]+b[1]*z3+b[2]*z3**2+b[3]*z3**3+b[4]*z3**4+b[5]*z3**5+b[6]*z3**6
    C1=(1+mibar*(8.*z3-2.*z3**2)/(1.-z3)**4+(1-mibar)*(20.*z3-27.*(z3**2)+12.*(z3**3)-2.*z3**4)/((1-z3)*(2.-z3))**2)**-1
    m2es3mat=np.outer(xi,xi)*np.outer(mi,mi)*(uij/T)*sigij**3
    m2e2s3mat=m2es3mat*(uij/T)
    m2es3,m2e2s3=np.sum(m2es3mat),np.sum(m2e2s3mat)
    fdisp=-2*np.pi*rho*I1*m2es3-np.pi*rho*mibar*C1*I2*m2e2s3
    
    # Association Contribution
    deltAiBj=gij*kapij*sigij**3*(np.exp(epsAiBj/T)-1.)
    rhoi=rho*xi
    def XAi_eq(XAi,rhoi,N,deltAiBj): return XAi-((1+np.sum(rhoi*XAi*N*deltAiBj.T,axis=1))**-1)
    XAiself=((-1.+np.sqrt(1.+4.*rho*np.diag(deltAiBj)))/(2.*rho*np.diag(deltAiBj)))
    XAiself[np.isnan(XAiself)]=1.
    XAi=wertheimiter(XAi_eq,XAiself.astype(xi.dtype),rhoi,N,deltAiBj)
    fassoc=np.sum((np.log(XAi)-1./2.*XAi+1./2.)*N*ntype*xi)#q=np.sum((np.log(XAi)-XAi+1.)*N*ntype*xi)-rho/2.*np.sum(np.outer(XAi,XAi)*deltAiBj*np.outer(N,N)*np.outer(xi,xi)*ntype)
    fres=fhc+fdisp+fassoc
    #_______________________dadx__________________________________
    # Hard Chain
    z3x=np.pi/6.*rho*mi*di**3
    z2x=np.pi/6.*rho*mi*di**2
    z1x=np.pi/6.*rho*mi*di**1
    z0x=np.pi/6.*rho*mi*di**0
    fhsx=(-z0x/z0*fhs+1./z0*(3.*(z1x*z2+z1*z2x)/(1.-z3)+3.*(z1*z2*z3x)/(1.-z3)**2+3.*(z2**2*z2x)/(z3*(1.-z3)**2)+(z2**3*z3x*(3.*z3-1.))/(z3**2*(1-z3)**3)+((3*z2**2*z2x*z3-2*z2**3*z3x)/z3**3-z0x)*np.log(1.-z3)+(z0-z2**3/z3**2)*z3x/(1.-z3)))
    gijx=(z3x/(1.-z3)**2+np.outer((dijmat),(3.*z2x/(1.-z3)**2+6.*z2*z3x/(1.-z3)**3))+np.outer((dijmat)**2,(4.*z2*z2x/(1.-z3)**3+6.*z2**2*z3x/(1.-z3)**4))).T.reshape((ncomp,ncomp,ncomp))
    giix=(z3x/(1.-z3)**2+np.outer((diimat),(3.*z2x/(1.-z3)**2+6.*z2*z3x/(1.-z3)**3))+np.outer((diimat)**2,(4.*z2*z2x/(1.-z3)**3+6.*z2**2*z3x/(1.-z3)**4)))
    fhcx=mi*fhs+mibar*fhsx-np.sum(xi*(mi-1)*gii**-1*giix.T,axis=1)-(mi-1)*np.log(gii) 
    
    # Dispersion Contribution
    m2es3x=2.*mi*np.sum(xi*mi*uij/T*sigij**3,axis=1)
    m2e2s3x=2.*mi*np.sum(xi*mi*(uij/T)**2*sigij**3,axis=1)
    C2=-C1**2*(mibar*(-4.*z3**2+20.*z3+8.)/(1.-z3)**5+(1-mibar)*(2.*z3**3+12.*z3**2-48.*z3+40.)/((1.-z3)*(2.-z3))**3)
    C1x=C2*z3x-C1**2*(mi*(8.*z3-2.*z3**2)/(1.-z3)**4-mi*(20.*z3-27.*(z3**2)+12.*(z3**3)-2.*z3**4)/((1.-z3)*(2.-z3))**2)
    ax=1./mibar**2*a0[1,:]+1./mibar**2*(3.-4./mibar)*a0[2,:]
    bx=1./mibar**2*b0[1,:]+1./mibar**2*(3.-4./mibar)*b0[2,:]
    I1x=(mi*(ax[0]+ax[1]*z3+ax[2]*z3**2+ax[3]*z3**3+ax[4]*z3**4+ax[5]*z3**5+ax[6]*z3**6)+a[1]*z3x+2.*a[2]*z3x*z3+3.*a[3]*z3x*z3**2+4.*a[4]*z3x*z3**3+5.*a[5]*z3x*z3**4+6.*a[6]*z3x*z3**5)
    I2x=(mi*(bx[0]+bx[1]*z3+bx[2]*z3**2+bx[3]*z3**3+bx[4]*z3**4+bx[5]*z3**5+bx[6]*z3**6)+b[1]*z3x+2.*b[2]*z3x*z3+3.*b[3]*z3x*z3**2+4.*b[4]*z3x*z3**3+5.*b[5]*z3x*z3**4+6.*b[6]*z3x*z3**5)
    fdispx=(-2.*np.pi*rho*(I1x*m2es3+I1*m2es3x)-np.pi*rho*((mi*C1*I2+mibar*C1x*I2+mibar*C1*I2x)*m2e2s3+mibar*C1*I2*m2e2s3x))
    
    # Association Contribution
    deltAiBjx=gijx*kapij*sigij**3*(np.exp(epsAiBj/T)-1.)
    qx=np.log(XAi)*N*ntype-rho/2*np.sum(np.sum(np.outer(XAi,XAi)*deltAiBjx*np.outer(N,N)*np.outer(xi,xi)*ntype,axis=2),axis=1)
    fresx=fhcx+fdispx+qx

    #_______________________dadrho__________________________________
    # Hard-Chain Contribution
    Zhs=z3/(1-z3)+3*z1*z2/(z0*(1-z3)**2)+(3*z2**3-z3*z2**3)/(z0*(1-z3)**3)
    gijz=z3/(1.-z3)**2+(dijmat)*(3.*z2/(1.-z3)**2+6.*z2*z3/(1.-z3)**3)+(dijmat)**2*(4.*z2*z2/(1.-z3)**3+6.*z2**2*z3/(1.-z3)**4)
    giiz=np.diag(gijz)
    Zhc=mibar*Zhs-np.sum(xi*(mi-1)*gii**-1*giiz)

    # Dispersion Contribution
    I1z=a[0]+2*a[1]*z3+3*a[2]*z3**2+4*a[3]*z3**3+5*a[4]*z3**4+6*a[5]*z3**5+7*a[6]*z3**6
    I2z=b[0]+2*b[1]*z3+3*b[2]*z3**2+4*b[3]*z3**3+5*b[4]*z3**4+6*b[5]*z3**5+7*b[6]*z3**6
    Zdisp=-2*np.pi*rho*(I1z*m2es3)-np.pi*rho*mibar*((C1*I2z+C2*z3*I2)*m2e2s3)

    # Association Contribution
    deltAiBjz=gijz*kapij*sigij**3*(np.exp(epsAiBj/T)-1.)
    Zassoc=-1./2.*np.sum((1.-XAi)*N*ntype*xi)-rho/2.*np.sum(np.outer(XAi,XAi)*deltAiBjz*np.outer(N,N)*np.outer(xi,xi)*ntype)
    Z=1+Zhc+Zdisp+Zassoc
    #_______________________mures__________________________________
    mures = fres + (Z-1.) + fresx - np.dot(xi, fresx)
    return fres,mures,Z
    

@njit(cache=True)
def etaiter(p,T,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB):
    def Z_obj(p,T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB):
        kB = 1.380649e-23
        di=sigi*(1.-0.12*np.exp(-3*ui/T))
        rho=6/np.pi*eta*(np.sum(mi*xi*di**3))**-1
        rhobar=rho*(10.**10)**3
        Zp=p/(rhobar*kB*T)
        _,_,Z1=ares(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB)
        return (Zp-Z1.real)
    def etaroots(fun,pp,p0,x,p1,p2,p3,p4,p5,p6,p7,p8,p9,tol=1E-8,iter=50):
        f=fun(pp,p0,x,p1,p2,p3,p4,p5,p6,p7,p8,p9)#.reshape(n)
        h = tol
        dx = h
        J= (fun(pp,p0,x+dx,p1,p2,p3,p4,p5,p6,p7,p8,p9)-f)/h
        for i in range(iter):
            if np.abs(f)<tol:
                return x 
            s=-1.*f/J
            x+=s
            df=fun(pp,p0,x,p1,p2,p3,p4,p5,p6,p7,p8,p9)-f
            J+=(df-J*s)/s
            f+=df
        return x
    eta0=0.45
    return etaroots(Z_obj,p,T,eta0,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB)

def vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N):
    etapures=[]
    for i in range(len(mi)):
        x=etaiter(p,T,np.asarray([1.]),np.asarray([mi[i]]),np.asarray([sigi[i]]),np.asarray([ui[i]]),np.asarray([epsAiBi[i]]),np.asarray([kapi[i]]),np.asarray([N[i]]),None,None)
        etapures.append(x)
    etapures=np.asarray(etapures)
    di=sigi*(1.-0.12*np.exp(-3*ui/T))
    NA = 6.0221407e23
    vmol=np.pi/6/etapures*mi*di**3/(10.**10)**3*NA
    return vmol

@njit(cache=True)
def SAFTSAC(T,vpure,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB):
    NA = 6.0221407e23
    #vpfracNET=(1-ksw*RH**2)/xi[0]
    #vmol=v0pNE/vpfracNET
    vmol=np.sum(vpure*xi)
    vpfrac=vpure/vmol
    di=sigi*(1.-0.12*np.exp(-3*ui/T))
    eta=np.pi/6*np.sum(mi*xi.real*di**3)/vmol/(10.**10)**3*NA
    etapure=np.pi/6*mi*di**3/vpure/(10.**10)**3*NA
    lngammaid=np.log(vpfrac)+1-vpfrac
    arespures=np.asarray([ares(T,val,np.asarray([1.]),np.asarray([mi[i]]),np.asarray([sigi[i]]),np.asarray([ui[i]]),np.asarray([epsAiBi[i]]),np.asarray([kapi[i]]),np.asarray([N[i]]),None,None)[0] for i,val in enumerate(etapure)])
    _,mures,Z1=ares(T,eta,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB)
    lngammares=mures-arespures
    lngammap=vpure/vmol*(Z1-1)
    return lngammaid+lngammares-lngammap

def THFaktor(T, vpure,xi,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB,idx=-1):
    nc = len(xi)
    h = 1E-26
    df = np.zeros([nc, nc])
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'complex128')
        dx[i] = h * 1j
        dx[idx] = - h * 1j  #x3+dx3=1-x1+dx1-x2+dx2=x3+dx1+dx2
        out =  SAFTSAC(T,vpure,xi+dx,mi,sigi,ui,epsAiBi,kapi,N,kij,kijAB)
        df[i] = out.imag/h
    return df.T   


#Test call, so all functions are compiled directly when it is imported
T=298.15
p=1E5
npoint=11
#Water
#XAiself =0.16 und 0.6
mi=np.asarray([1.20469,2.38269789])
sigi=np.asarray([2.797059952,3.1771])
ui=np.asarray([353.95,198.24])
epsAiBi=np.asarray([2425.67,2653.4])
kapi=np.asarray([0.04509,0.032384])
N=np.asarray([1.,1.])
x1=np.linspace(0,1,npoint)
x2=1-x1
xi=np.vstack((x1,x2))
vpures=vpure(p,T,mi,sigi,ui,epsAiBi,kapi,N)
lngammai=np.asarray([SAFTSAC(T,vpures,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N,None,None).flatten() for i,val in enumerate(xi[0,:])])
Gammai=np.asarray([THFaktor(T,vpures,np.ascontiguousarray(xi[:,i]),mi,sigi,ui,epsAiBi,kapi,N,None,None).flatten() for i,val in enumerate(xi[0,:])])