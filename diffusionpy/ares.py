import numpy as np
from numba import njit

@njit(['Tuple((f8, f8[:], f8))(f8,f8,f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1])',
        'Tuple((c16, c16[:], c16))(f8,c16,c16[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1])'],cache=True)
def ares(T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA):
    """calculate the reduced residual helmholtz energy, the chemical potential and the real gas factor

    Args:
        T (float): temperature
        xi (array_like): mole/mass fraction. Becomes the mass fraction when the molar mass Mi is not None
        mi (array_like): segment number
        si (array_like): segment diameter
        ui (array_like): dispersion energy
        eAi (array_like): association energy
        kAi (array_like): association volume
        NAi (array_like): association sites (only symmetric)
        kij (array_like): Matrix of binary interaction parameters for dispersion.
        kijA (array_like): Matrix of binary interaction parameters for association.
    Return:
        ares (array_like): reduced residual helmholtz energy
        mures (array_like): reduced residual chemical potential
        Zres (aaray_like): real gas factor
    """
    def np_add_outer(a):  
        """create a outer product but with addition instead of multiplication""" 
        return a.reshape(len(a),1)+a
    def wertheim_iter(fun,x,p1,p2,p3,p4):
        """solve this non-linear equation system for the nonbonded association sites XAi""" 
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
    def mat(vec):
        m=np.zeros((ncomp,ncomp))
        for i,j,k in zip(*np.triu_indices(ncomp,k=1),range(len(vec)+1)): m[i,j]=vec[k];m[j,i]=vec[k]
        return m
    #Initializeki
    # for i,j,k in zip(*np.triu_indices(ncomp,k=1),range(ncomp)): kijmat[i,j]=kij[k]
    # kijAmat+=kijAmat.T-np.diag(kijAmat)
    #Mixing rules
    mibar=np.sum(xi*mi)
    di=si*(1.-0.12*np.exp(-3*ui/T))
    rho=6/np.pi*eta*(np.sum(mi*xi.real*di**3))**-1 
    dij=np_add_outer(di)
    sij=np_add_outer(si)/2.

    kijmat=mat(kij)
    kijAmat=mat(kijA)
    uij=np.outer(ui,ui)**0.5*(1.-kijmat)
    kAij=np.outer(kAi,kAi)**0.5*(np.outer(si,si)**0.5/sij)**3
    eAij=np_add_outer(eAi)/2.*(1.-kijAmat)

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
    def XAi_eq(XAi,xi,rho,NAi,deltAij):
        """solve this non-linear equation system for the nonbonded association sites XAi"""
        return XAi-((1+np.sum(xi*rho*XAi*NAi*deltAij.T,axis=1))**-1)
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
    mures = fres + (Z-1) + fresx - np.dot(xi, fresx)
    return fres,mures,Z
    
