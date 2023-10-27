import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from numba import njit,config,prange
import time
from .FEM_collocation import collocation,collocation_space

# config.DISABLE_JIT = True





def Diffusion_MS(tint,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,dlnai_dlnwi=None,**kwargs):
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
    Examples:
        >>> t=np.linspace(0,300,100)
        >>> Dvec=np.asarray([1E-13])
        >>> wi0=np.asarray([0.3,0.7])
        >>> wi8=np.asarray([0.7,0.3])
        >>> L=1E-6
        >>> Mi=np.asarray([18.,72.])
        >>> mobile=np.asarray([True,False])
        >>> wt=Diffusion_MS(tint,L,Dvec,wi0,wi8,Mi,mobile)
    See Also:
        diffusionpy.D_Matrix
    """
    @njit(['f8[::1](f8, f8[:],f8[:], f8[:,:,::1,::1], i8[::1], i8[::1], f8[::1], f8[:,:], b1, f8[::1],f8[:,:],f8[:,::1])'],cache=True)
    def ode(t,x,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB):
        """change in the weight fraction with time"""
        def vanishing_check(dwvdt,wv):# Realitätscheck gegen neg. Konz.
            s1,s2=wv.shape
            for i in range(s1):
                for j in range(s2):
                    if (dwvdt[i,j]<0 and wv[i,j]<0 ): dwvdt[i,j]=0
                    if (dwvdt[i,j]>0 and wv[i,j]>1 ): dwvdt[i,j]=0  
                return dwvdt
        def np_linalg_solve(A,b):
            """solve a Batch of linear system of equations"""
            ret = np.zeros_like(b)
            for i in range(b.shape[1]):
                ret[:, i] = np.linalg.solve(A[:,:,i],b[:,i])
            return ret
        refsegment=np.argmin(Mi)
        ri= Mi/Mi[refsegment]
        def BIJ_Matrix(D,wi,mobiles):
            """create the friction matrix which needs to be invertable"""
            nc,nz_1=wi.shape
            B=np.zeros((nc,nc,nz_1))
            for i in range(nc):
                Dii=np.zeros(nz_1)
                for j in range(nc):
                    if j!=i:
                        Dij=-ri[j]*wi[i,:]/D[i,j] 
                        B[i,j,:]=Dij
                        Dii+=ri[i]*wi[j,:]/D[i,j]
                B[i,i,:]=Dii
            return B[mobiles,:,:][:,mobiles,:]
        nTH,nz_1=dmuext.shape
        wv=np.reshape(np.ascontiguousarray(x),(nTH,nz_1))    
        nc=len(Mi)
        wi=np.zeros((nc,nz_1))
        wi[mobiles,:]=wv
        wi0_immobiles=np.expand_dims(wi0[immobiles]/np.sum(wi0[immobiles]),axis=1)
        wi[immobiles,:]=(1-np.expand_dims(np.sum(wv,axis=0),axis=0))*wi0_immobiles
        for j in range(nc):
            wi[j,-1]=np.interp(t,tint,wiB[:,j])
        wv[:,-1]=wi[mobiles,-1]
        dwv=np.zeros((nTH,nz_1))
        for j in range(nTH):
            dwv[j,:]=collocation(wv[j,:],nz_1,True)
        dwv[:,0]=0
        dlnwv=dwv/wv
        THFaktor_=np.zeros((nz_1,nTH,nTH))
        for i in range(nz_1):
            for j in range(nTH):
                for k in range(nTH):
                    THFaktor_[i,j,k]=np.interp(t,tint,THFaktor[:,i,j,k])
        dlnav=np.zeros_like(dlnwv)
        for i in range(nz_1): 
            dlnav[:,i]=THFaktor_[i,...]@np.ascontiguousarray(dlnwv[:,i]) #if np.linalg.det(THFaktor_[i,...])>0.001 else dlnwv[:,i]#-(THFaktor_[i,...]@np.ascontiguousarray(dlnwv[:,i]))
        B=BIJ_Matrix(D,wi,mobiles)
        
        if allflux:
            for i in range(nz_1):
                B[:,:,i]=B[:,:,i]+1/np.max(D)*np.max(ri)*np.outer(wv[:,i],np.ones_like(wv[:,i]))
        dmuv=dlnav+dmuext
        omega=(np.sum(wi0[immobiles],axis=0)/(1-np.sum(wv,axis=0)))**-1 if not allflux else np.ones(nz_1)
        dv=wv*dmuv*omega #if swelling else wv*dmuv/np.atleast_2d(ri).T*omega              
        jv=np_linalg_solve(B,dv) #if not allflux else np_linalg_solve(B,dv[:-1,:])
        djv=np.zeros((nTH,nz_1))
        for j in range(nTH):
            djv[j,:]=collocation(jv[j,:],nz_1,False)    
        dwvdt=np.zeros_like(djv)
        for j in range(nTH):
            dwvdt[j,:]=djv[j,:]-np.sum(djv,axis=0)*wv[j,:] if not allflux else djv[j,:]
        return  (dwvdt*omega).flatten()#vanishing_check(djv,rhov).flatten()


    print("------------- Initialization and postprocessing ----------------")
    start1=time.time_ns()
    nc=len(wi0)
    nz=kwargs["nz"] if "nz" in kwargs else 20
    z,nE=collocation_space(nz+1)
    dz=L/(nE)
    D=D_Matrix(Dvec/dz**2,nc)    
    zvec=z*L
    nf=int(np.sum(mobile))
    nt=len(tint)
    allflux=nc==nf
    nTH=nf #if not allflux else nc-1
    mobiles=np.where(mobile)[0] #if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] #if not allflux else np.asarray([nc-1],dtype=np.int64)
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    wiinit=(wi0*np.ones((nz+1,nc))).T
    wiinit[:,-1]=wi8
    wvinit=wiinit[mobiles,:]
    #Construct TH Factor
    THFaktor=np.asarray([[np.eye(nTH)]*(nz+1)]*nt)
    if dlnai_dlnwi is not None:
        if len(dlnai_dlnwi.shape)==2:
            dlnai_dlnwi=dlnai_dlnwi[None,:,:]*np.ones((nt,nc,nc))
        if len(dlnai_dlnwi.shape)==3:
            dlnai_dlnwi=dlnai_dlnwi[:,None,:,:]*np.ones((nt,nz+1,nc,nc))
        if len(dlnai_dlnwi.shape)==4:
            slc1=np.ix_(np.asarray(range(nt)),np.asarray(range(nz+1)),mobiles, mobiles) 
            slc2=np.ix_(np.asarray(range(nt)),np.asarray(range(nz+1)),immobiles, mobiles)
            massbalancecorrection=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles[None,None,:,None],axis=2) 
            THFaktor=dlnai_dlnwi[slc1]-massbalancecorrection[:,:,None,:]
            # THFaktorave=np.average(np.average(THFaktor,axis=0),axis=0)
            # THFaktor=THFaktorave[None,None,:,:]*np.ones((nt,nz+1,nTH,nTH))
    xinit=wvinit.flatten()
    dmuext=np.zeros((nTH,nz+1))
    wiB=kwargs['witB'] if "witB" in kwargs else wi8[None,:]*np.ones((nt,nc))
    drhovdtB=np.zeros(nTH)
    return_stress=False
    return_alpha=False
            
    if "EJ" in kwargs or "etaJ" in kwargs or "exponent" in kwargs: 
        from .relaxation import relaxation_mode
        rho0i=kwargs["rho0i"] if "rho0i" in kwargs else 1200.*np.ones(nc)
        xinit,ode=relaxation_mode(wvinit,ode,kwargs["EJ"],kwargs["etaJ"],kwargs["exponent"],Mi[mobiles],1/rho0i[mobiles],tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB)
        return_stress=True

    if "deltHSL" in kwargs or "TSL" in kwargs or "cpSL" in kwargs  or "DAPI" in kwargs  or "sigma" in kwargs  or "kt" in kwargs or "g" in kwargs or "crystallize" in kwargs: 
        from .crystallization import crystallization_mode
        lngi_tz=interp1d(tint,kwargs["lngi_tz"],axis=0,bounds_error=False)
        rho0i=kwargs["rho0i"] if "rho0i" in kwargs else 1200.*np.ones(nc)
        xinit,ode=crystallization_mode(wvinit,ode,mobiles,immobiles,kwargs["crystallize"],wi0,wi8,rho0i,Mi,kwargs["deltaHSL"],kwargs["TSL"],kwargs["cpSL"],kwargs["DAPI"],kwargs["sigma"],kwargs["kt"],kwargs["g"],lngi_tz)
        return_alpha=True
    if "kappaii" in kwargs:
        from .liquidseperation import liquidseperation_mode
        kappaii=kwargs["kappaii"]/dz**2
        xinit,ode=liquidseperation_mode(wvinit,ode,kappaii,tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB)

    print("------------- Start diffusion modeling ----------------")
    start=time.time_ns()
    # sol=nbkode.BDF5(ode,tint[0],xinit,params=[tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB]).run(tint)
    # sol=nbrk_ode(ode,(tint[0],tint[-1]),xinit,args=(tint,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,rhoiB,drhovdtB),rk_method=0,t_eval=tint)#rtol=1E-2,atol=1E-3)
    sol=solve_ivp(ode,(tint[0],tint[-1]),xinit,args=(tint,THFaktor,mobiles,immobiles,Mi,D,allflux,wi0,dmuext,wiB),method="Radau",t_eval=tint,atol=1E-3)#rtol=1E-2,atol=1E-3)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    if not sol["success"]: raise Exception(sol["message"]+f" The time step of failing was {tint[len(sol['y'])]} seconds")# the initial conditions are returned instead ----------------"); #return wi0*np.ones((nc,nt)).T 
    # x_sol=np.exp(sol["y"] ) if "EJ" not in kwargs else sol["y"] 
    x_sol=sol["y"]
    # rhoend=x_sol[:(nz+1)*nTH,-1]
    
    nt=len(tint)
    wt=np.zeros((nc,nt))
    wik=np.zeros((nt,nc,nz+1))

    for k in range(nt):
        wvk=np.reshape(x_sol[:(nz+1)*nTH,k],(nTH,nz+1))
        wik[k,:,:]= wiinit
        wik[k,mobiles,:]=wvk
        wik[k,immobiles,:]=(1-np.sum(wik[k,mobiles,:],axis=0))*wi0_immobiles[:,None]
        wik[k,:,-1]=wiB[k,:]
        wt[:,k]=np.sum(wik[k,:,:-1]/nz,axis=1)
    Lt=np.sum(wi0[immobiles])/(1-np.sum(wt[mobiles,:],axis=0))*L
    # plt.plot(t,rhoik[:,1,-1])
    end1=time.time_ns()
    print(f"------------- Initialization and postprocessing took {((end1-start1)-(end-start))/1E9} seconds----------------")
    
    if return_stress:
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
    nd=(nc-1)*nc//2 #Anzahl der Diffusionskoeffizienten
    if len(Dvec)!=nd: #Länge des Vektors der Diff.Koeffs. muss mit Anzahl der Koeffs- übereinstimmen
        raise Exception("Wrong length of ``Dvec``. Provide array with ``(nc-1)*nc/2`` entries") #Sonst wird ein Fehler ausgegeben
    else:
        D=np.zeros((nc,nc)) # 0-Matrix der Dimension nc,nc wird erstellt
        D[np.triu_indices_from(D,k=1)]=Dvec #Dreiecksmatrix mit Werten aus Dvec wird erstellt
        D=D.T # D wird transformiert
        D[np.triu_indices_from(D,k=1)]=Dvec #Dreiecksmatrix mit Werten aus Dvec wird erstellt
    return D # D wird zurückgegeben

def wegstein(fun,x):
    """Solving via wegsteins method"""
    tol=1E-6 #Toleranzbereich
    maxiter=10 #Anzahl max. Iterationen
    f = fun(x)
    xx=f    
    dx = xx - x
    ff = fun(xx)
    df=ff-f
    for i in range(maxiter): #Iteration von 0 bis maxiter
        e=np.linalg.norm(dx.flatten(),2)/np.prod(x.shape)
        print(f"iter {i+1}: ||F|| = {e}")
        if e<tol: 
            return xx 
        a = df/dx
        # q= np.nan_to_num(a/(a-1),nan=0,posinf=0,neginf=0)
        q=0 # resort to fixed point iteration as it works best for this case
        x=xx
        xx = q * xx + (1-q) * ff
        # xx[:,-1]=1-np.sum(xx[:,:-1],axis=1)
        f=ff
        ff = fun(xx)    
        df=ff-f
        dx=xx-x
    return xx  

def Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,dlnai_dlnwi_fun=None,**kwargs):
    """iterates dlnai_dlnwi as a function of the concentration wi
    See also:
        diffusionpy.Diffusion_MS
    
    """
    _,wt_old,_,_=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,**kwargs,full_output=True)
    _,_,nz_1=wt_old.shape
    def wt_fix(wtz):
        
        
        print("------------- Start PC-SAFT modeling ----------------")
        wtz=np.ascontiguousarray(np.swapaxes(wtz,1,2))
        start=time.time_ns()
        dlnai_dlnwi=dlnai_dlnwi_fun(wtz)
        end=time.time_ns()
        print("------------- PC-SAFT modeling took "+str((end-start)/1E9)+" seconds ----------------")
        return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,full_output=True,**kwargs)[1]
    wtopt=wegstein(wt_fix,wt_old)
    wtopt=np.ascontiguousarray(np.swapaxes(wtopt,1,2))
    dlnai_dlnwiopt=dlnai_dlnwi_fun(wtopt)
    return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=full_output,dlnai_dlnwi=dlnai_dlnwiopt,**kwargs)


