import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root,fixed_point
from scipy.interpolate import interp1d
# from scipy.special import logsumexp
from numba import njit,config,prange
import time

# config.DISABLE_JIT = True

@njit(['f8[::1](f8, f8[:,::1],f8[:], f8[:,:,::1,::1], i8[::1], i8[::1], f8[::1], f8[:,:], b1, b1, f8, f8[::1],f8[:,:],f8[:,::1],f8[::1])'],cache=True)
def drhodt(t,rhov,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB):
    """change in the partial density with time"""
    def vanishing_check(drhovdt,rhov):
        s1,s2=rhov.shape
        for i in range(s1):
            for j in range(s2):
                if (drhovdt[i,j]<0 and rhov[i,j]<1E-4 ): drhovdt[i,j]=0 
            return drhovdt
    def averaging(a):
        """make rolling average over a 2D array"""
        return (a[1:,...]+a[:-1,...])/2.
    
    def np_linalg_solve(A,b):
        """solve a Batch of linear system of equations"""
        ret = np.zeros_like(b)
        for i in range(b.shape[1]):
            ret[:, i] = np.linalg.solve(A[:,:,i], b[:,i])
        return ret
    def BIJ_Matrix(D,wi,mobiles,allflux):
        """create the friction matrix which needs to be invertable"""
        nc,nz_1=wi.shape
        B=np.zeros((nc,nc,nz_1))
        for i in range(nc):
            Din=wi[-1,:]/D[i,-1] if (i+1)!=nc else np.zeros(nz_1)
            Dii=np.zeros(nz_1)
            for j in range(nc):
                if j!=i:
                    Dij=-wi[i,:]/D[i,j] 
                    if allflux: Dij+=Din
                    B[i,j,:]=Dij
                    Dii+=wi[j,:]/D[i,j]
            if allflux: Dii+=Din
            B[i,i,:]=Dii
        return B[mobiles,:,:][:,mobiles,:] if not allflux else B[:-1,:-1,:]


    refsegment=np.argmin(Mi)
    ri= Mi[mobiles]/Mi[refsegment]
    nc=len(Mi)
    nTH,nz_1=rhov.shape
    #_________
    nP=2

    nE=(nz_1-1)//(nP-1) # 5 Collocation points used to calculate the amount of finite elements there should be an error when these do not match the given nz
    # S=np.asarray([
    #     [-25./3., 16., -12., 16./3., -1.],
    #     [-1., -10./3., 6., -2., 1./3.],
    #     [1/3, -8./3., 0., 8./3, -1./3.], 
    #     [-1./3., 2., -6., 10./3., 1.],
    #     [1., -16./3., 12., -16., 25./3.]])
    # S=np.asarray([[-9.183300132670368, 13.95870750793668, -7.32835571279449, 3.5529483375281643, -0.9999999999999911],  # Legendre collocation matrix calculated for 5 collocation points
    # [-1.7403293111129658, -1.3744088030005464, 4.3546484316145255, -1.682881139058434, 0.44297082155741346], 
    # [0.5458250331675965, -2.601439792602136, 3.552713678800501e-15, 2.6014397926021253, -0.5458250331675947], 
    # [-0.4429708215574084, 1.6828811390584317, -4.354648431614518, 1.3744088030005417, 1.7403293111129672], 
    # [0.9999999999999858, -3.552948337528136, 7.328355712794462, -13.95870750793668, 9.183300132670368]])*nE

    #chebyshev 1st order
    S=np.asarray([[-9.47213595499959, 13.708203932499401, -6.472135954999608, 3.2360679774998005, -1.],
    [-2, -0.8541019662496925, 4., -1.6180339887499053, 0.4721359549995813], 
    [0.6180339887498896, -2.618033988749879, 0., 2.6180339887499073, -0.6180339887498976],
    [-0.4721359549995922, 1.618033988749914, -4., 0.8541019662496961, 2.], 
    [1., -3.2360679774998857, 6.472135954999686, -13.70820393249943, 9.472135954999601]])
    #4th order
    # S=np.asarray([[-5.828427124746187, 8.242640687119275, -3.414213562373087, 0.9999999999999973], [-1.4142135623730971, -0.41421356237309154, 2.414213562373092, -0.5857864376269044], [0.5857864376269034, -2.4142135623730896, 0.414213562373088, 1.4142135623730976], [-0.9999999999999956, 3.414213562373087, -8.242640687119275, 5.828427124746187]])
    S=np.asarray([[-1.0, 1.0], [-1.0, 1.0]]) #basically finite differences
    # S=np.asarray([[-19.19566935808956, 28.293504037135108, -14.591793886480787, 8.987918414870705, -5.603875471610536, 3.10991626417524, -1.0000000000001705],
    # [-3.603875471609483, -1.960771339502083, 8.09783467904494, -4.000000000000474, 2.317667207394258, -1.2469796037176133, 0.3961245283903642],
    # [0.8900837358252135, -3.8780021506949454, -0.5211062828031919, 4.98791841487004, -2.246979603717616, 1.1099162641748106, -0.34183037765438584],
    # [-0.44504186791264555, 1.5549581320873287, -4.048917339522259, -3.552713678800501e-15, 4.048917339522426, -1.5549581320874104, 0.44504186791263667],
    # [0.3418303776543146, -1.1099162641749691, 2.246979603717414, -4.987918414869792, 0.5211062828031414, 3.878002150694952, -0.8900837358251975], 
    # [-0.3961245283905829, 1.246979603718002, -2.3176672073945834, 4.000000000000232, -8.097834679045107, 1.9607713395018953, 3.603875471609564],
    # [1.0000000000002558, -3.109916264176036, 5.603875471611104, -8.987918414870933, 14.591793886479763, -28.293504037134312, 19.19566935808919]])
    # S=np.asarray([[-3.0, 4.0, -1.0],
    # [-1.0, 0.0, 1.0], 
    # [1.0, -4.0, 3.0]])
        #_________
    # if t<tint[0] or t>tint[-1]: print(f"Im outside of the range wit {t} seconds")
    rhoi=np.zeros((nc,nz_1))
    rhoi[mobiles,:]=rhov
    rhoi[immobiles,:]=rho*np.expand_dims(wi0[immobiles],1) if not allflux else np.expand_dims(rho-np.sum(rhoi[:-1,:],axis=0),0)
    for j in range(nc):
        rhoi[j,-1]=np.interp(t,tint,rhoiB[:,j])
    
    rhov[:,-1]=rhoi[mobiles,-1]
    wi=rhoi/np.sum(rhoi,axis=0)
    wv=wi[mobiles,:]


    #_________
    dwv=np.zeros((nTH,nz_1))
    
    for i in range(nE):
        P0=i*(nP-1)
        P8=(i+1)*(nP-1)+1
        wvP=wv[:,P0:P8]
        dwvP=np.zeros_like(wvP)
        if i==0: dwvP_=np.zeros_like(wvP)
        for j in range(nTH):
            dwvP[j,:]=(S@np.ascontiguousarray(wvP[j,:]))
        if i>0: dwvP[:,0]=dwvP_[:,-1]
        dwv[:,P0:P8]=dwvP
        dwvP_=dwvP
        
    dwv[:,0]=0
    dlnwv=dwv/wv
    # dlnwv=drhov*(1-wv)/rhov
    # dlnwv[:,0]=0
    #_________
    

    THFaktor_=np.zeros((nz_1,nTH,nTH))
    for i in range(nz_1):
        for j in range(nTH):
            for k in range(nTH):
                THFaktor_[i,j,k]=np.interp(t,tint,THFaktor[:,i,j,k])
    # if THFaktor.ndim==2:
    #     dlnai=dlnwv+np.diff(THFaktor[:,mobiles],axis=0).T
    # else:
    dlnav=np.zeros_like(dlnwv)
    for i in range(nz_1): 
        dlnav[:,i]=THFaktor_[i,...]@np.ascontiguousarray(dlnwv[:,i])
    B=BIJ_Matrix(D,wi,mobiles,allflux)
    dmuv=dlnav+dmuext
    omega=rho/np.sum(rhoi,axis=0)
    # di=rho*wvbar*dmui/np.atleast_2d(ri).T*omega if swelling else rhovbar*dmui/np.atleast_2d(ri).T
    dv=rho*wv*dmuv/np.atleast_2d(ri).T*omega if swelling else rhov*dmuv/np.atleast_2d(ri).T
    jv=np_linalg_solve(B,dv) #if not allflux else np_linalg_solve(B,dv[:-1,:])
    # if allflux:
    #     nonetflux=np.expand_dims(-np.sum(jv,axis=0),0)
    #     jv=np.vstack((jv,nonetflux))
    
    # jv[:,0]=0

    djv=np.zeros((nTH,nz_1))
    for i in range(nE-1,-1,-1):
        P0=i*(nP-1)
        P8=(i+1)*(nP-1)+1
        jvP=jv[:,P0:P8]
        djvP=np.zeros_like(jvP)
        if i==(nE-1): djvP_=np.zeros_like(jvP)
        for j in range(nTH):
            djvP[j,:]=(S@np.ascontiguousarray(jvP[j,:]))
        if i<(nE-1): djvP[:,-1]=djvP_[:,0]   
        djv[:,P0:P8]=djvP
        djvP_=djvP
    #_________
    # dji=np.diff(np.hstack((jiB,ji)))
    djv[:,-1]=0

    # drhovdt=np.hstack((dji,np.atleast_2d(drhovdtB).T))
    return  vanishing_check(djv,rhov).flatten()


def Diffusion_MS(tint,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,dlnai_dlnwi=None,swelling=False,**kwargs):
    """
    Method that computes the multi-component diffusion kinetics 
    
    Args:
        t (array_like): time
        L (float) : dry thickness /m
        Dvec (array_like): Vector of diffusion coefficients. See diffusionpy.D_Matrix                       /m^2/s
        wi0 (array_like): Mass fractions at t=0               /-
        wi8 (array_like): Mass fraction at t=infinity         /-
        Mi (array_like):   Molar mass of components nc         /g/mol
        mobile(array_like): boolean vector indicating the mobility of a component
        dlnai_dlnwi (array_like): estimate for DlnaiDlnx at t          /-
        Keyword Arguments:
            wiB (array_like): Hello \n
            rho0iB (array_like): Hello \n
    Returns:
        ndarray:   
        if ``full_output=False``: \n
        Matrix ``wt`` of mass fractions at t /- \n

        if ``full_output=True``: \n 
        Matrix of mass fractions at t       /- \n
        Matrix of mass fractions at t,z     /- \n
    See Also:
        diffusionpy.D_Matrix
    """
    print("------------- Initialization and postprocessing ----------------")
    start1=time.time_ns()
    nc=len(wi0)
    
    nz=kwargs["nz"] if "nz" in kwargs else 20
    nP=2
    nE=(nz)//(nP-1)
    dz=L/(nE)
    D=D_Matrix(Dvec/dz**2,nc)    
    # NS=np.asarray([0.,0.20289048,0.5,0.79710952,1.]) #Legendre
    # NS=np.asarray([0.,0.25,0.5,0.75,1.]) #equidistant
    NS=np.asarray([0.,0.19098301,0.5,0.80901699,1.]) #chebyshev 1st order
    # NS=np.asarray([0.0, 0.2928932188134524, 0.7071067811865476, 1.0]) #chebyshev 1st order 4th polynomial
    NS=np.asarray([0.,1.])
    # NS=np.asarray([0.0, 0.09903113209758088, 0.27747906604368555, 0.5, 0.7225209339563144, 0.9009688679024191, 1.0])
    # NS=np.asarray([0.0, 0.5, 1.0])
    z=np.asarray([0.])
    zE=np.linspace(0,1,nE+1)
    for i in range(nE):
        z0=zE[i]
        z8=zE[i+1]
        zP=NS*(z8-z0)+z0
        z=np.hstack((z,zP[1:]))
    zvec=z*L
    nf=int(np.sum(mobile))
    nt=len(tint)
    rho=1200. if "rho0i" not in kwargs else np.sum(kwargs["rho0i"]*wi0)
    if "rho0i" not in kwargs  : rho0i=rho*np.ones(nc)
    if "rho0i" in kwargs : rho0i=kwargs["rho0i"]
    allflux=nc==nf
    nTH=nf if not allflux else nc-1
    mobiles=np.where(mobile)[0] if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] if not allflux else np.asarray([-1],dtype=np.int64)
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    rhoiinit=(rho*wi0*np.ones((nz+1,nc))).T
    rhovinit=rhoiinit[mobiles,:]
    #Construct TH Factor
    THFaktor=np.asarray([[np.eye(nTH)]*(nz+1)]*nt)
    if dlnai_dlnwi is not None:
        if len(dlnai_dlnwi.shape)==3:
            # slc1=np.ix_(np.asarray(range(nt)),mobiles, mobiles) #if not allflux else np.ix_(np.asarray(range(nt)),np.arange(0,nc-1,dtype=np.int64), np.arange(0,nc-1,dtype=np.int64)) 
            # slc2=np.ix_(np.asarray(range(nt)),immobiles, mobiles) #if not allflux else np.ix_(np.asarray(range(nt)),np.asarray([nc-1]), np.arange(0,nc-1,dtype=np.int64)) 
            # massbalancecorrection=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles[None,:,None],axis=1) #if not allflux else np.sum(dlnai_dlnwi[slc2],axis=1)
            # THFaktor_=dlnai_dlnwi[slc1]-massbalancecorrection[:,None,:]
            # THFaktor_= interp1d(tint,THFaktor_,axis=0,bounds_error=False,fill_value=(THFaktor_[0,:,:],THFaktor_[-1,:,:]))
            # THFaktor= lambda t: np.ones((nz,nTH,nTH))*THFaktor_(t)
            dlnai_dlnwi=np.swapaxes(np.asarray([dlnai_dlnwi]*(nz+1)),0,1)
            
        if len(dlnai_dlnwi.shape)==4:
            slc1=np.ix_(np.asarray(range(nt)),np.asarray(range(nz+1)),mobiles, mobiles) #if not allflux else np.ix_(np.asarray(range(nt)),np.asarray(range(nz+1)),np.arange(0,nc-1,dtype=np.int64), np.arange(0,nc-1,dtype=np.int64)) 
            slc2=np.ix_(np.asarray(range(nt)),np.asarray(range(nz+1)),immobiles, mobiles) #if not allflux else np.ix_(np.asarray(range(nt)),np.asarray(range(nz+1)),np.asarray([nc-1]), np.arange(0,nc-1,dtype=np.int64)) 
            massbalancecorrection=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles[None,None,:,None],axis=2) #if not allflux else np.sum(dlnai_dlnwi[slc2],axis=2)
            THFaktor=dlnai_dlnwi[slc1]-massbalancecorrection[:,:,None,:]
            # THFaktor= interp1d(t,THFaktor_,axis=0,bounds_error=False,fill_value=(THFaktor_[0,:,:,:],THFaktor_[-1,:,:,:]))
            # THFaktor= lambda t: np.ones((nz,nTH,nTH))*THFaktor_(t)
    #Set up the PDEs boundary,initial conditions and additional driving forces.
    #____________________________________
    #def default_mode():
    xinit=rhovinit.flatten()
    dmuext=np.zeros((nTH,nz+1))
    rhoiB=wi8[None,:]*rho*np.ones((nt,nc))
    drhovdtB=np.zeros(nTH)
    return_sigma=False
    return_alpha=False
    # float64, array(float64, 2d, C), array(float64, 1d, C), array(float64, 4d, C), array(int64, 1d, C), array(int64, 1d, C), array(float64, 1d, C), array(float64, 2d, F), bool, bool, float64, array(float64, 1d, C), array(float64, 2d, C), array(float64, 2d, C), array(float64, 1d, C)"
    @njit(['f8[::1](f8, f8[:],f8[:], f8[:,:,::1,::1], i8[::1], i8[::1], f8[::1], f8[:,:], b1, b1, f8, f8[::1],f8[:,:],f8[:,::1],f8[::1])'],cache=True)
    def ode(t,x,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB):
        """create generic ode function for drhodt"""
        rhov=np.reshape(np.ascontiguousarray(x),(nTH,nz+1))
        return drhodt(t,rhov,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB)
        
    # Mechanical equation of state (MEOS)      
    if "EJ" in kwargs or "etaJ" in kwargs or "exponent" in kwargs: 
        from .relaxation import relaxation_mode
        xinit,ode=relaxation_mode(rhovinit,ode,kwargs["EJ"],kwargs["etaJ"],kwargs["exponent"],Mi[mobiles],1/rho0i[mobiles],tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB)
        return_sigma=True

    if "deltHSL" in kwargs or "TSL" in kwargs or "cpSL" in kwargs  or "DAPI" in kwargs  or "sigma" in kwargs  or "kt" in kwargs or "g" in kwargs or "crystallize" in kwargs: 
        from .crystallization import crystallization_mode
        lngi_tz=interp1d(tint,kwargs["lngi_tz"],axis=0,bounds_error=False)
        xinit,ode=crystallization_mode(rhovinit,ode,mobiles,immobiles,kwargs["crystallize"],wi0,wi8,rho0i,Mi,kwargs["deltaHSL"],kwargs["TSL"],kwargs["cpSL"],kwargs["DAPI"],kwargs["sigma"],kwargs["kt"],kwargs["g"],lngi_tz)
        return_alpha=True
        # THFaktor=lngi_tz
        # THFaktor= lambda t: (np.ones((nz,nTH,nTH))-alpha)*THFaktor_(t)
    #_____________________________________
    # if "lngi_tz" in kwargs:
    #     lngi_tz=interp1d(t,kwargs["lngi_tz"],axis=0,bounds_error=False,kind="quadratic")
    #     THFaktor=lngi_tz

    if "witB" in kwargs: rhoiB=kwargs['witB']*rho#interp1d(tint,kwargs['witB']*rho,axis=0,bounds_error=False,fill_value=(kwargs['witB'][0]*rho,kwargs['witB'][-1]*rho))
    def wrapode(t,x,ode,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB):
        """evaluate time dependent functions and insert into the generic odes"""
        THFaktor=THFaktor(t)    if callable(THFaktor)   else THFaktor
        #hier könnte das eins minus alpha irgendwo muss das x ausgepackt werden
        rhoiB=rhoiB(t)          if callable(rhoiB)      else rhoiB
        drhovdtB=drhovdtB(t,x)  if callable(drhovdtB)   else drhovdtB
        dxdt=ode(t,x,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB)
        return dxdt.flatten()
    print("------------- Start diffusion modeling ----------------")
    start=time.time_ns()
    # sol=nbkode.BDF5(ode,tint[0],xinit,params=[tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB]).run(tint)
    # sol=nbrk_ode(ode,(tint[0],tint[-1]),xinit,args=(tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB),rk_method=0,t_eval=tint)#rtol=1E-2,atol=1E-3)
    sol=solve_ivp(ode,(tint[0],tint[-1]),xinit,args=(tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB),method="Radau",t_eval=tint)#rtol=1E-2,atol=1E-3)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    if not sol["success"]: raise Exception(sol["message"])# the initial conditions are returned instead ----------------"); #return wi0*np.ones((nc,nt)).T 
    # x_sol=np.exp(sol["y"] ) if "EJ" not in kwargs else sol["y"] 
    x_sol=sol["y"]
    # rhoend=x_sol[:(nz+1)*nTH,-1]
    
    nt=len(tint)
    wt=np.zeros((nc,nt))
    wik=np.zeros((nt,nc,nz+1))
    rhoik=np.zeros((nt,nc,nz+1))
    rhok=np.zeros(nt)
    rhoit=np.zeros((nt,nc))
    for k in range(nt):
        rhovk=np.reshape(x_sol[:(nz+1)*nTH,k],(nTH,nz+1))
        rhoik[k,:,:]=rhoiinit
        
        rhoik[k,mobiles,:]=rhovk
        rhoik[k,:,-1]=rhoiB[k,:]
        # if allflux: rhoik[k,immobiles,:]=rho-np.sum(rhoik[k,:-1,:],axis=0) 
        for i in range(nc):
            wik[k,i,:]=rhoik[k,i,:]/np.sum(rhoik[k,:,:],axis=0)
        wik[k,immobiles,-1]=(1-np.sum(wik[k,mobiles,-1],axis=0))*wi0_immobiles
        rhok[k]=np.sum(np.sum(rhoik[k,:,:-1]/nz,axis=1),axis=0)
        rhoit[k,:]=np.sum(rhoik[k,:,:-1]/nz,axis=1)
        wt[:,k]=np.sum(wik[k,:,:-1]/nz,axis=1)
    wt=(rhoit/rhok[:,None]).T
    Lt=rhok/rho*L
    # plt.plot(t,rhoik[:,1,-1])
    end1=time.time_ns()
    print(f"------------- Initialization and postprocessing took {((end1-start1)-(end-start))/1E9} seconds----------------")
    
    if return_sigma:
        nJ=len(kwargs["EJ"])
        sigmaJ=np.reshape(x_sol[(nz+1)*nTH:(nz+1)*(nTH+nJ),:],(nz+1,nJ,nt))
        return(wt.T,sigmaJ) if not full_output else (wt.T,wik,zvec,Lt,sigmaJ)

    elif return_alpha:
        alpha=x_sol[(nz+1)*nTH:(nz+1)*(nTH+1),:]
        r=x_sol[(nz+1)*(nTH+1):(nz+1)*(nTH+2),:]
        return (wt.T,alpha,r) if not full_output else (wt.T,wik,zvec,Lt,alpha,r)

    else:
        return wt.T if not full_output else (wt.T,wik,zvec,Lt)


def D_Matrix(Dvec,nc):
    """
    Creates a symmetric square Matrix ``Dmat`` of dimension ``nc`` 
    using the elements in the vector ``Dvec``.
    It is assumed that the elements of ``Dvec`` fit precisly
    in the upper and lower triangular entries of ``Dmat``.
    Args:
        Dvec (array_like): Must have the length of  ``(nc-1)*nc/2`` to fit in the diagionals of the result matrix
        nc (int): Dimension of ``Dmat``.
    Returns:
        ndarray: square matrix ``Dmat`` of shape ``(nc,nc)``
    Raises:
        Exception
            Wrong length of ``Dvec``. Provide array with ``(nc-1)*nc/2`` entries
    Examples:
        >>> Dvec=np.array([1E-13,2E-13,3E-13])
        >>> nc=3
        >>> Dmat=D_Matrix(Dvec,nc)
        >>> Dmat
        array([[0.e+00, 1.e-13, 2.e-13],
               [1.e-13, 0.e+00, 3.e-13],
               [2.e-13, 3.e-13, 0.e+00]])
    """
    nd=(nc-1)*nc//2 
    if len(Dvec)!=nd: 
        raise Exception("Wrong length of ``Dvec``. Provide array with ``(nc-1)*nc/2`` entries")
    else:
        D=np.zeros((nc,nc))
        D[np.triu_indices_from(D,k=1)]=Dvec
        D=D.T
        D[np.triu_indices_from(D,k=1)]=Dvec
    return D

def wegstein(fun,x):
    """Solving via wegsteins method"""
    tol=1E-6
    maxiter=50
    f = fun(x)
    xx=f    
    dx = xx - x
    ff = fun(xx)
    df=ff-f
    for i in range(maxiter):
        e=np.linalg.norm(dx.flatten(),2)/np.prod(x.shape)
        print(f"iter {i+1}: ||F|| = {e}")
        if e<tol: 
            return xx 
        a = df/dx
        q= np.nan_to_num(a/(a-1),nan=0,posinf=0,neginf=0)
        x=xx
        xx = q * xx + (1-q) * ff
        # xx[:,-1]=1-np.sum(xx[:,:-1],axis=1)
        f=ff
        ff = fun(xx)    
        df=ff-f
        dx=xx-x
    return xx  

def Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,swelling=False,dlnai_dlnwi_fun=None,**kwargs):
    """iterates dlnai_dlnwi as a function of the concentration wi
    See also:
        diffusionpy.Diffusion_MS
    
    """
    nt=len(t)
    nc=len(wi0)
    
    # if dlnai_dlnwi_fun is None : Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output,None,swelling,**kwargs)
    # dlnai_dlnwi=np.stack([dlnai_dlnwi_fun(wi8*0.5+wi0*0.5)]*nt)
    # wt_old=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,**kwargs)
    _,wt_old,_,_=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=swelling,**kwargs,full_output=True)
    _,_,nz_1=wt_old.shape
    def wt_obj(wt_old):
        wtz=wt_old.reshape((nt,nc,nz_1))
        # wtz=(wtz[:,:,:1]+wtz[:,:,:-1])/2
        dlnai_dlnwi=np.asarray([[dlnai_dlnwi_fun(col) for col in row.T] for row in wtz])
        residual=(wt-Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,**kwargs)).flatten()/nt
        return residual
    def wt_fix(wtz):
        # wtz=(wtz[:,:,:1]+wtz[:,:,:-1])/2
        wtz=np.ascontiguousarray(np.swapaxes(wtz,1,2))
        print("------------- Start PC-SAFT modeling ----------------")
        start=time.time_ns()
        # dlnai_dlnwi=np.asarray([[dlnai_dlnwi_fun(np.ascontiguousarray(col)) for col in row.T] for row in wtz])
        dlnai_dlnwi=dlnai_dlnwi_fun(wtz)
        end=time.time_ns()
        print("------------- PC-SAFT modeling took "+str((end-start)/1E9)+" seconds ----------------")
        return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,full_output=True,**kwargs)[1]
    method=kwargs["method"] if "method" in kwargs  else "df-sane"
    if method=="df-sane":
        wtopt=root(wt_obj,wt_old.flatten(),method=method,options={"disp":True,"maxfev":30,"fatol":1E-6})["x"].reshape((nt,nc)) 
    elif method=="wegstein":
        wtopt=wegstein(wt_fix,wt_old)
    else:
        for i in range(kwargs["maxit"]):
            wtopt=wt_fix(wt_old)
            wt_old=wtopt
    if full_output:
        # dlnai_dlnwiopt=np.asarray([[dlnai_dlnwi_fun(col) for col in row.T] for row in (wtopt[:,:,:1]+wtopt[:,:,:-1])/2])
        # wtopt=(wtopt[:,:,:1]+wtopt[:,:,:-1])/2
        wtopt=np.ascontiguousarray(np.swapaxes(wtopt,1,2))
        dlnai_dlnwiopt=dlnai_dlnwi_fun(wtopt)
        return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=full_output,dlnai_dlnwi=dlnai_dlnwiopt,swelling=swelling,**kwargs)
    else:
        return np.average(wtopt,axis=2)

def convert(x,M,axis=0):
    """convert fractions. e.g mass fractions into mole fractions"""
    y=np.empty(x.shape,M.dtype)
    np.copyto(y.T,M)
    return x*y/np.sum(x*y,axis=axis)

def correctMB(dlnai_dlnwi,mobile,nt,wi0):
    mobiles=np.where(mobile)[0]
    immobiles=np.where(~mobile)[0]
    allflux=len(mobile)==np.sum(mobile)
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    wi0_mobiles=wi0[mobiles]/np.sum(wi0[mobiles])
    slc1=np.ix_(np.asarray(range(nt)),mobiles, mobiles) if not allflux else (np.asarray(range(nt)),np.arange(0,nc-1,dtype=np.int64), np.arange(0,nc-1,dtype=np.int64)) 
    slc2=np.ix_(np.asarray(range(nt)),immobiles, mobiles) if not allflux else (np.asarray(range(nt)),-1, np.arange(0,nc-1,dtype=np.int64)) 
    slc3=np.ix_(np.asarray(range(nt)),mobiles, immobiles) if not allflux else (np.asarray(range(nt)),-1, np.arange(0,nc-1,dtype=np.int64)) 
    massbalancecorrection1=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles[None,:,None],axis=1) if not allflux else np.sum(dlnai_dlnwi[slc2],axis=1)
    massbalancecorrection2=np.sum(dlnai_dlnwi[slc3]*wi0_mobiles[None,:,None],axis=1) if not allflux else np.sum(dlnai_dlnwi[slc3],axis=1)
    THFaktor_=dlnai_dlnwi[slc1]-np.stack((massbalancecorrection1[:,None,:],massbalancecorrection2[:,None,:]),axis=1)
    return THFaktor_