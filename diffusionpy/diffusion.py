import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import root
from numba import njit,config,prange
import time
from .PCSAFT import dlnai_dlnxi_loop,dlnai_dlnxi,SAFTSAC,vpure
from .surface import time_dep_surface


def Diffusion_MS(tint,L,Dvec,wi0,wi8,mobile,T=298.15,p=1E5,saftpar=None,**kwargs):
    """
    Method that computes the multi-component diffusion kinetics 
    
    Args:
        t (array_like): time
        L (float) : dry thickness /m
        Dvec (array_like): Vector of diffusion coefficients. See diffusionpy.D_Matrix                       /m^2/s
        wi0 (array_like): Mass fractions at t=0               /-
        wi8 (array_like): Mass fraction at t=infinity         /-
        mobile(array_like): boolean vector indicating the mobility of a component
        dlnai_dlnwi (array_like): estimate for DlnaiDlnx at t          /-
        Keyword Arguments:
            wiB (array_like): weight fraction at z=L \n
    Returns:
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
        >>> wt=Diffusion_MS(tint,L,Dvec,wi0,wi8,mobile)
    See Also:
        diffusionpy.D_Matrix
    """
    # @njit(['f8[::1](f8, f8[:],f8[:], f8[:,:,::1,::1], i8[::1], i8[::1], f8[:,:], b1, f8[:,:],f8[:,:],f8[:,::1])'],cache=True)
    def ode(t,x,tint,THFaktor,mobiles,immobiles,D,allflux,wi0,dmuext,wiB):
        def BIJ_Matrix(D,wi,mobiles):
            """create the friction matrix containning the diffusion coefficients"""
            nc,nz_1=wi.shape
            B=np.zeros((nc,nc,nz_1))
            for i in range(nc):
                Dii=np.zeros(nz_1)
                for j in range(nc):
                    if j!=i:
                        Dij=-wi[i,:]/D[i,j] 
                        B[i,j,:]=Dij
                        Dii+=wi[j,:]/D[i,j]
                B[i,i,:]=Dii
            return B[mobiles,:,:][:,mobiles,:]
        nTH,nz_1=dmuext.shape
        wv=np.reshape(np.ascontiguousarray(x),(nTH,nz_1))    
        nc,_=wi0.shape
        wi=np.zeros((nc,nz_1))
        wi[mobiles,:]=wv
        wi0_immobiles=wi0[immobiles,:]/np.sum(wi0[immobiles,:],axis=0)
        wi[immobiles,:]=(1-np.expand_dims(np.sum(wv,axis=0),axis=0))*wi0_immobiles
        for j in range(nc): wi[j,-1]=np.interp(t,tint,wiB[:,j])
        wv[:,-1]=wi[mobiles,-1]
        dwv=np.zeros((nTH,nz_1))
        dwi=np.zeros((nc,nz_1))
        for j in range(nTH): dwv[j,1:]=np.diff(wv[j,:])
        dwv[:,0]=0
        dwi[mobiles,:]=dwv
        dwi[immobiles,:]=-np.sum(dwv,axis=0)*wi0_immobiles
        THFaktor_=np.zeros((nz_1,nc,nc))
        for i in range(nz_1):
            for j in range(nc):
                for k in range(nc):
                    THFaktor_[i,j,k]=np.interp(t,tint,THFaktor[:,i,j,k])
        dai=np.zeros_like(dwi)
        for i in range(nz_1): dai[:,i]=THFaktor_[i,:,:]@dwi[:,i]
        dav=dai[mobiles,:]
        B=BIJ_Matrix(D,wi,mobiles)
        if allflux:
            for i in range(nz_1):
                B[:,:,i]=B[:,:,i]+1/np.max(D)*np.outer(wv[:,i],np.ones_like(wv[:,i]))
        dmuv=dav+dmuext*wv
        omega=(np.sum(wi0[immobiles,:],axis=0)/(1-np.sum(wv,axis=0)))**-1 if not allflux else np.ones(nz_1)
        dv=dmuv*omega
        jv = np.zeros_like(dv)
        for i in range(jv.shape[1]):  jv[:, i] = np.linalg.solve(B[:,:,i],dv[:,i]) 
        jv[:,0]=0.
        djv=np.zeros((nTH,nz_1))
        for j in range(nTH):
            djv[j,:-1]=np.diff(jv[j,:])   
        djv[:,-1]=0.
        dwvdt=np.zeros_like(djv)
        for j in range(nTH): dwvdt[j,:]=djv[j,:]-np.sum(djv,axis=0)*wv[j,:] if not allflux else djv[j,:]
        return  (dwvdt*omega).flatten()#vanishing_check(djv,rhov).flatten()


    start1=time.time_ns()
    nc=len(wi0)
    nz=kwargs["nz"] if "nz" in kwargs else 20
    z=np.linspace(0,L,nz+1)
    dz=L/(nz)
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
    


    xinit=wvinit.flatten()
    dmuext=np.zeros((nTH,nz+1))
    wiB=time_dep_surface(tint,wi0,wi8,mobile,kwargs['taui']) if "taui" in kwargs else wi8[None,:]*np.ones((nt,nc))
    drhovdtB=np.zeros(nTH)
    return_stress=False
    return_alpha=False

    if "EJ" in kwargs or "etaJ" in kwargs or "exponent" in kwargs: 
        from .relaxation import relaxation_mode
        ode=relaxation_mode(ode,**kwargs)
        nJ=len(kwargs["EJ"])
        sigmaJ0=np.zeros((nz+1,nJ))
        sigmaJB=np.zeros((nJ))
        sigmaJ0[-1,:]=sigmaJB
        xinit=np.hstack((wvinit.flatten(),sigmaJ0.flatten()))
        return_stress=True
    if "A" in kwargs or "B" in kwargs or "C" in kwargs  or "n" in kwargs: 
        from .crystallization import crystallization_mode
        alpha0=np.zeros(nz+1)
        # r0=1E-20*np.ones(nz+1)
        # xinit=np.hstack((wvinit.flatten(),alpha0.flatten(),r0.flatten()))
        xinit=np.hstack((wvinit.flatten(),alpha0.flatten()))
        # kwargs["lnS"]=interp1d(tint,kwargs["lnS"],axis=0,bounds_error=False)
        wv_fun=kwargs['wv_fun'] if 'wv_fun' in kwargs else None
        ode=crystallization_mode(ode,kwargs['A'],kwargs['B'],kwargs['n'],T,saftpar,wv_fun)
        return_alpha=True
    
    def getsol(sol):
        if not sol["success"]: raise Exception(sol["message"]+f" The time step of failing was {tint[len(sol['y'])]} seconds")# the initial conditions are returned instead ----------------"); #return wi0*np.ones((nc,nt)).T  
        x_sol=sol["y"]
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
        if return_stress:
            nJ=len(kwargs["EJ"])
            sigmaJ=np.reshape(x_sol[(nz+1)*nTH:(nz+1)*(nTH+nJ),:],(nz+1,nJ,nt))
            return (wt.T,wik,zvec,Lt,sigmaJ)
        elif return_alpha:
            alpha=x_sol[(nz+1)*nTH:(nz+1)*(nTH+1),:]
            return (wt.T,wik,zvec,Lt,alpha)
        else:    
            return (wt.T,wik,zvec,Lt)
    
    def solveodes(wik=None):
        THFaktor=np.asarray([[np.eye(nc)]*(nz+1)]*nt) 
        if saftpar is not None:
            if wik is not None:
                wik=np.ascontiguousarray(np.swapaxes(wik,1,2))
                if 'vpure' not in saftpar: saftpar['vpure']=vpure(p,T,**saftpar)
                THFaktor=np.nan_to_num(Gammaij(T,wik,saftpar))
        sol=solve_ivp(ode,(tint[0],tint[-1]),xinit,args=(tint,THFaktor,mobiles,immobiles,D,allflux,wi0[:,None]*np.ones((nc,nz+1)),dmuext,wiB),method="Radau",t_eval=tint)#rtol=1E-2,atol=1E-3)


        # plt.plot(t,rhoik[:,1,-1])

        return getsol(sol)

    solopt=wegstein(solveodes,solveodes()[1]) if (saftpar is not None) and ('A' not in kwargs)  else solveodes()
    return solopt


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


def massbalancecorrection(dlnai_dlnwi,wi0,mobile):
    points=dlnai_dlnwi.shape[:-2] if len(dlnai_dlnwi.shape)>2 else []
    mobiles=np.where(mobile)[0] 
    immobiles=np.where(~mobile)[0]
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    slc1=np.ix_(*([np.arange(val) for val in  points]+[mobiles, mobiles]))
    slc2=np.ix_(*([np.arange(val) for val in  points]+[mobiles, immobiles]))
    correction=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles[...,None,:],axis=-1)  
    THFaktor=(dlnai_dlnwi[slc1]-correction[...,:,None])#*wi[...,None,mobiles]
    return THFaktor

def Gammaij(T,wi,par):
    Mi=par["Mi"]
    ri=Mi/np.min(Mi)
    dlnai_dlnwi=dlnai_dlnxi_loop(T,np.ascontiguousarray(wi),**par) if len(wi.shape)==3 else dlnai_dlnxi(T,np.ascontiguousarray(wi),**par)
    with np.errstate(divide='ignore',invalid='print'):
        dlnai_dlnwi=dlnai_dlnwi/wi[...,None,:]*wi[...,:,None]/ri[...,:,None]
    return dlnai_dlnwi


def DIdeal2DReal(Dvec,wave,wi0,dlnai_dlnwi,mobile,realtoideal=False):
    nc=wi0.shape[0]
    nf=int(np.sum(mobile))
    allflux=nc==nf
    nTH=nf if not allflux else nc-1
    mobiles=np.where(mobile)[0] 
    THFaktor=massbalancecorrection(dlnai_dlnwi,wi0,mobile)
    if allflux: THFaktor+=np.max(np.linalg.eig(THFaktor)[0])*np.outer(wave,np.ones_like(wave))   
    if realtoideal: THFaktor=np.linalg.inv(THFaktor)
    def BIJ(D,wi,mobiles):
        nc=wi.shape[0]
        B=np.zeros((nc,nc))
        for i in range(nc):
            Dii=0
            for j in range(nc):
                if j!=i:
                    Dij=-wi[i]/D[i,j] 
                    B[i,j]=Dij
                    Dii+=wi[j]/D[i,j]
            B[i,i]=Dii
        if allflux: B+=1/np.max(D)*np.outer(wi,np.ones_like(wi))   
        return B[mobiles,:][:,mobiles]
    D=D_Matrix(Dvec,nc)
    C=BIJ(D,wave,mobiles)

    # B=(THFaktor/wave[None,mobiles])@C #somehow in codrying almost identical
    B=THFaktor@C #symmetric without correction
    # B=THFaktor@C
    def Bopt(Dsoll):
        D=D_Matrix(Dsoll,nc)
        xist=BIJ(D,wave,mobiles).flatten()
        xsoll=B.flatten()
        return np.sum((1-xist/xsoll)**2)
    from scipy.optimize import minimize
    opt=minimize(Bopt,Dvec,method="Nelder-Mead",bounds=(((-1E-1,1E-1),)*len(Dvec)))

    DMS=opt["x"]
    DMS[DMS<0]=1E-5
    return DMS



def wegstein(fun,x):
    """Solving via wegsteins method"""
    tol=1E-6 #Toleranzbereich
    maxiter=20 #Anzahl max. Iterationen
    sol = fun(x)
    f=sol[1]
    xx=f    
    dx = xx-x
    ssol = fun(xx)
    ff=ssol[1]
    df=ff-f
    for i in range(maxiter): #Iteration von 0 bis maxiter
        e=np.linalg.norm(df.flatten(),2)/np.prod(x.shape)
        # print(f"iter {i+1}: ||F|| = {e}")
        if e<tol: 
            return ssol
        with np.errstate(divide='ignore',invalid='print'):    
            a = df/dx
            q= np.average(np.fmin(np.fmax(np.nan_to_num(a/(a-1),nan=0,posinf=0,neginf=0),0),1))
        # q=0
        # print(f"iter {i+1}: q = {q}")
        x=xx
        xx = q * xx + (1-q) * ff
        f=ff
        ssol = fun(xx)
        ff=ssol[1]    
        df=ff-f
        dx = xx - x
    return ssol  

