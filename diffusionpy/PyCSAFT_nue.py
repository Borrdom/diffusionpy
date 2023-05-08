import numpy as np
from numba import njit,config
# config.DISABLE_JIT=True
# import os
# os.environ['NUMBA_DEBUG_CACHE'] = "1"

# def logging_jit(func):
    
#     def inner(*args, **kwargs):
#         origsigs = set(func.signatures)
#         result = func(*args, **kwargs)
#         newsigs = set(func.signatures)
#         if newsigs != origsigs:
#             new = (newsigs ^ origsigs).pop()
#              # PUT YOUR LOGGER HERE!
#             print("Numba compiled for signature: {}".format(new))
#         return result
#     return inner

# @logging_jit
@njit(['Tuple((f8, f8[:], f8))(f8,f8,f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])',
        'Tuple((c16, c16[:], c16))(f8,c16,c16[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])'],cache=True)
def ares(T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA):
    def np_add_outer(a): 
        return a.reshape(len(a),1)+a
    def wertheim_iter(fun,x,p1,p2,p3,p4):
        tol=1E-8
        iter=50
        n=len(x)
        f=fun(x,p1,p2,p3,p4)
        h = tol
        J=np.ones((n,n)).astype(p1.dtype)
        for i in range(n):
            dx = x*0
            dx[i] = h
            J[:,i] = (fun(x + dx,p1,p2,p3,p4)-f)/h
        for i in range(iter):
            if np.linalg.norm(f,2)<tol:
                return x 
            s=np.linalg.solve(J+np.eye(n)*tol,-1.*f)
            ff=fun(x+s,p1,p2,p3,p4)
            lamb=1.
            for j in range(iter):
                if np.linalg.norm(ff,2)*(1.-lamb*1E-1)<np.linalg.norm(f,2): break 
                lamb*=1./4.
                s*=lamb 
                ff=fun(x+s,p1,p2,p3,p4) 
            x+=s
            df=fun(x,p1,p2,p3,p4)-f
            J+=np.outer((df-np.dot(J,s)),s)/np.dot(s,s)
            f+=df
        return x
    npoly_B=7
    npoly_A=3
    a0=np.ones((npoly_A,npoly_B))
    b0=np.ones((npoly_A,npoly_B))
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


    #Initializeki
    # for i,j,k in zip(*np.triu_indices(ncomp,k=1),range(ncomp)): kijmat[i,j]=kij[k]
    # kijAmat+=kijAmat.T-np.diag(kijAmat)
    #Mixing rules
    mibar=np.sum(xi*mi)
    di=si*(1.-0.12*np.exp(-3*ui/T))
    rho=6/np.pi*eta*(np.sum(mi*xi.real*di**3))**-1 
    dij=np_add_outer(di)
    sij=np_add_outer(si)/2.
    uij=np.outer(ui,ui)**0.5*(1.-kij)
    kAij=np.outer(kAi,kAi)**0.5*(np.outer(si,si)**0.5/sij)**3
    eAij=np_add_outer(eAi)/2.*(1.-kijA)

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
    m2es3mat=np.outer(xi,xi)*np.outer(mi,mi)*(uij/T)*sij**3
    m2e2s3mat=m2es3mat*(uij/T)
    m2es3,m2e2s3=np.sum(m2es3mat),np.sum(m2e2s3mat)
    fdisp=-2*np.pi*rho*I1*m2es3-np.pi*rho*mibar*C1*I2*m2e2s3
    
    # Association Contribution
    deltAij=gij*kAij*sij**3*(np.exp(eAij/T)-1.)
    rhoi=rho*xi
    def XAi_eq(XAi,xi,rho,NAi,deltAij): return XAi-((1+np.sum(xi*rho*XAi*NAi*deltAij.T,axis=1))**-1)
    deltAi=np.fmax(np.diag(deltAij),1E-300)
    XAi0=((-NAi+np.sqrt(NAi**2+4.*NAi*rho*deltAi))/(2.*rho*deltAi))
   #XAi0[np.isnan(XAi0)]=1.
    XAi=wertheim_iter(XAi_eq,XAi0.astype(xi.dtype),xi,rho,NAi,deltAij)
    fassoc=np.sum((np.log(XAi)-1./2.*XAi+1./2.)*NAi*ntype*xi)#q=np.sum((np.log(XAi)-XAi+1.)*NAi*ntype*xi)-rho/2.*np.sum(np.outer(XAi,XAi)*deltAij*np.outer(NAi,NAi)*np.outer(xi,xi)*ntype)
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
    m2es3x=2.*mi*np.sum(xi*mi*uij/T*sij**3,axis=1)
    m2e2s3x=2.*mi*np.sum(xi*mi*(uij/T)**2*sij**3,axis=1)
    C2=-C1**2*(mibar*(-4.*z3**2+20.*z3+8.)/(1.-z3)**5+(1-mibar)*(2.*z3**3+12.*z3**2-48.*z3+40.)/((1.-z3)*(2.-z3))**3)
    C1x=C2*z3x-C1**2*(mi*(8.*z3-2.*z3**2)/(1.-z3)**4-mi*(20.*z3-27.*(z3**2)+12.*(z3**3)-2.*z3**4)/((1.-z3)*(2.-z3))**2)
    ax=1./mibar**2*a0[1,:]+1./mibar**2*(3.-4./mibar)*a0[2,:]
    bx=1./mibar**2*b0[1,:]+1./mibar**2*(3.-4./mibar)*b0[2,:]
    I1x=(mi*(ax[0]+ax[1]*z3+ax[2]*z3**2+ax[3]*z3**3+ax[4]*z3**4+ax[5]*z3**5+ax[6]*z3**6)+a[1]*z3x+2.*a[2]*z3x*z3+3.*a[3]*z3x*z3**2+4.*a[4]*z3x*z3**3+5.*a[5]*z3x*z3**4+6.*a[6]*z3x*z3**5)
    I2x=(mi*(bx[0]+bx[1]*z3+bx[2]*z3**2+bx[3]*z3**3+bx[4]*z3**4+bx[5]*z3**5+bx[6]*z3**6)+b[1]*z3x+2.*b[2]*z3x*z3+3.*b[3]*z3x*z3**2+4.*b[4]*z3x*z3**3+5.*b[5]*z3x*z3**4+6.*b[6]*z3x*z3**5)
    fdispx=(-2.*np.pi*rho*(I1x*m2es3+I1*m2es3x)-np.pi*rho*((mi*C1*I2+mibar*C1x*I2+mibar*C1*I2x)*m2e2s3+mibar*C1*I2*m2e2s3x))
    
    # Association Contribution
    deltAijx=gijx*kAij*sij**3*(np.exp(eAij/T)-1.)
    qx=np.log(XAi)*NAi*ntype-rho/2*np.sum(np.sum(np.outer(XAi,XAi)*deltAijx*np.outer(NAi,NAi)*np.outer(xi,xi)*ntype,axis=2),axis=1)
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
    deltAijz=gijz*kAij*sij**3*(np.exp(eAij/T)-1.)
    Zassoc=-1./2.*np.sum((1.-XAi)*NAi*ntype*xi)-rho/2.*np.sum(np.outer(XAi,XAi)*deltAijz*np.outer(NAi,NAi)*np.outer(xi,xi)*ntype)
    Z=1+Zhc+Zdisp+Zassoc
    #_______________________mures__________________________________
    mures = fres + (Z-1.) + fresx - np.dot(xi, fresx)
    return fres,mures,Z
    

#@njit(cache=True)
def eta_iter(p,T,xi,mi,si,ui,eAi,kAi,NAi,kij=np.asarray([[0.]]),kijA=np.asarray([[0.]])):
    def Z_obj(p,T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij=np.asarray([[0.]]),kijA=np.asarray([[0.]])):
        kB = 1.380649e-23
        di=si*(1.-0.12*np.exp(-3*ui/T))
        rho=6/np.pi*eta*(np.sum(mi*xi*di**3))**-1
        rhobar=rho*(10.**10)**3
        Zp=p/(rhobar*kB*T)
        _,_,Z1=ares(T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA)
        return (Zp-Z1.real)
    def eta_roots(fun,pp,p0,x,p1,p2,p3,p4,p5,p6,p7,p8,p9,tol=1E-8,iter=50):
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
    return eta_roots(Z_obj,p,T,eta0,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA)

def vpure(p,T,mi,si,ui,eAi,kAi,NAi,**kwargs):
    etapures=[]
    for i in range(len(mi)):
        x=eta_iter(p,T,np.asarray([1.]),np.asarray([mi[i]]),np.asarray([si[i]]),np.asarray([ui[i]]),np.asarray([eAi[i]]),np.asarray([kAi[i]]),np.asarray([NAi[i]]),np.asarray([[0.]]),np.asarray([[0.]]))
        etapures.append(x)
    etapures=np.asarray(etapures)
    di=si*(1.-0.12*np.exp(-3*ui/T))
    NA = 6.0221407e23
    vmol=np.pi/6/etapures*mi*di**3/(10.**10)**3*NA
    return vmol

#@njit(cache=True)
def lngi(T,vpure,xi,mi,si,ui,eAi,kAi,NAi,Mi=None,kij=np.asarray([[0.]]),kijA=np.asarray([[0.]]),**kwargs):
    NA = 6.0221407e23
    xi=np.ascontiguousarray(xi)
    #vpfracNET=(1-ksw*RH**2)/xi[0]
    #vmol=v0pNE/vpfracNET
    wi=xi
    if Mi is not None: xi=xi/Mi/np.sum(xi/Mi)
    vmol=np.sum(vpure*xi)
    vpfrac=vpure/vmol
    di=si*(1.-0.12*np.exp(-3*ui/T))
    eta=np.pi/6*np.sum(mi*xi.real*di**3)/vmol/(10.**10)**3*NA
    etapure=np.pi/6*mi*di**3/vpure/(10.**10)**3*NA
    lngi_id=np.log(vpfrac)+1-vpfrac
    arespures=np.asarray([ares(T,val,np.asarray([1.]),np.asarray([mi[i]]),np.asarray([si[i]]),np.asarray([ui[i]]),np.asarray([eAi[i]]),np.asarray([kAi[i]]),np.asarray([NAi[i]]),np.asarray([[0.]]),np.asarray([[0.]]))[0] for i,val in enumerate(etapure)])
    _,mures,Z1=ares(T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA)
    lngi_res=mures-arespures
    lngi_p=vpure/vmol*(Z1-1)
    return lngi_id+lngi_res-lngi_p+np.log(xi/wi)

#@njit(cache=True)
def lnphi_TP(p,T,xi,mi,si,ui,eAi,kAi,NAi,Mi=None,kij=np.asarray([[0.]]),kijA=np.asarray([[0.]]),**kwargs):
    xi=np.ascontiguousarray(xi)
    if Mi is not None: xi=xi/Mi/np.sum(xi/Mi)
    etamix=np.asarray([eta_iter(p,T,np.ascontiguousarray(xi[:,i]),mi,si,ui,eAi,kAi,NAi,kij,kijA) for i,val in enumerate(xi[0,:])])
    lnphi=np.asarray([ares(T,val,np.ascontiguousarray(xi[:,i]),mi,si,ui,eAi,kAi,NAi,kij,kijA)[1].flatten() for i,val in enumerate(etamix)])
    return lnphi

def dlnai_dlnxi(T, vpure,xi,mi,si,ui,eAi,kAi,NAi,Mi=None,kij=np.asarray([[0.]]),kijA=np.asarray([[0.]]),idx=None,**kwargs):
    xi=np.ascontiguousarray(xi)
    nc = len(xi)
    h = 1E-26
    df = np.zeros([nc, nc])
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'c16')
        dx[i] = h * 1j
        if idx is not None: dx[idx] = - h * 1j #x3+dx3=1-x1+dx1-x2+dx2=x3+dx1+dx2
        #wi_= wi+dx
        #xi_=wi_/Mi/(wi_/Mi).sum()
        out =  lngi(T,vpure,xi+dx,mi,si,ui,eAi,kAi,NAi,Mi,kij,kijA)
        df[i] = out.imag/h
    return df*xi+np.eye(nc)


T=298.15
p=1E5
npoint=2
#Water
#XAi0 =0.16 und 0.6

par={"mi":np.asarray([1.20469,2.38269789]),
 "si":np.asarray([2.797059952,3.1771]),
 "ui":np.asarray([353.95,198.24]),
 "eAi":np.asarray([2425.67,2653.4]),
 "kAi":np.asarray([0.04509,0.032384]),
 "NAi":np.asarray([1.,1.])}
x1=np.linspace(0,1,npoint)
x2=1-x1
xi=np.vstack((x1,x2))
vpures=vpure(p,T,**par)
lngiammai=np.asarray([lngi(T,vpures,xi[:,i],**par).flatten() for i,val in enumerate(xi[0,:])])
Gammai=np.asarray([dlnai_dlnxi(T,vpures,xi[:,i],**par).flatten() for i,val in enumerate(xi[0,:])])