import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.interpolate import interp1d#,InterpolatedUnivariateSpline
from numba import njit,config,stencil
import time
config.DISABLE_JIT = True


@njit(['f8[:,:](f8, f8[:,::1], f8[:,::1], i8[::1], i8[::1], f8[::1], f8[:,:], b1, b1, f8, f8[::1],f8[:,:],f8[::1],f8[::1])'],cache=True)
def drhodt(t,wv,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,wiB,dwvdtB):
    """change in the partial density with time"""

    def np_linalg_solve(A,b):
        """solve a Batch of linear system of equations with numba"""
        ret = np.empty_like(b)
        for i in range(b.shape[1]):
            ret[:, i] = np.linalg.solve(A[:,:,i], b[:,i])
        return ret

    def np_gradient(f,fourth_order=False):
        """reimplementation of numpy gradient with second order over the entire domain for use with numba.  Optional third order for the center"""
        out = np.empty_like(f, np.float64)
        out[1:-1,:] = (f[2:,:] - f[:-2,:]) / 2.0
        out[0,:] = -(3*f[0,:] - 4*f[1,:] + f[2,:]) / 2.0
        out[-1,:] = (3*f[-1] - 4*f[-2] + f[-3]) / 2.0

        if fourth_order:
            out[2:-2,:] *=4./3.
            out[2:-2,:] -= (f[3:,:] - f[:-3,:]) / 12.0
        return out
    
    def BIJ_Matrix(D,wi,mobiles,allflux):
        """create the friction matrix which needs to be invertable"""
        nc,nz=wi.shape
        B=np.zeros((nc,nc,nz))
        for i in range(nc):
            Din=wi[-1,:]/D[i,-1] if (i+1)!=nc else np.zeros(nz)
            Dii=np.zeros(nz)
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
    nTH,nz_1=wv.shape
    wi=np.zeros((nc,nz_1))
    wi[mobiles,:]=wv
    wi[immobiles,:]=np.expand_dims(wi0[immobiles],1)
    wi[immobiles,-1]=wiB[immobiles]
    if not np.any(dwvdtB): wi[mobiles,-1]=wiB[mobiles]

    dlnai=THFaktor@(np_gradient(np.log(wi[mobiles,:]).T).T)
    B=BIJ_Matrix(D,wi,mobiles,allflux) 
    dmui=dlnai+dmuext
    di=wi[mobiles,:]*dmui/np.atleast_2d(ri).T
    ji=np_linalg_solve(B,di)
    jiB=np.zeros(nTH)
    ji[:,0]=jiB
    dji=np_gradient(ji.T).T
    dji[:,-1]=dwvdtB
    dwvdt=(dji-np.sum(dji,0)*wi[mobiles,:])
    # dwvdt=dji
    return  dwvdt


def Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,dlnai_dlnwi=None,swelling=False,**kwargs):
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
    nc=len(wi0)
    
    nz=kwargs["nz"] if "nz" in kwargs else 20
    dz=L/nz
    D=D_Matrix(Dvec/dz**2,nc)
    zvec=np.linspace(0,L,nz+1)
    nf=np.sum(mobile)
    nt=len(t)
    rho=1200. if "rho0i" not in kwargs else np.sum(kwargs["rho0i"]*wi0)
    if "rho0i" not in kwargs  : rho0i=rho*np.ones(nc)
    if "rho0i" in kwargs : rho0i=kwargs["rho0i"]
    allflux=nc==nf
    nTH=nf if not allflux else nc-1
    mobiles=np.where(mobile)[0] if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] if not allflux else np.asarray([-1],dtype=np.int64)
    wi0_immobiles=wi0[immobiles]/np.sum(wi0[immobiles])
    wiinit=(wi0*np.ones((nz+1,nc))).T
    wvinit=wiinit[mobiles,:]
    #Construct TH Factor
    THFaktor=np.eye(nTH)
    if dlnai_dlnwi is not None:
        slc1=np.ix_(np.asarray(range(nt)),mobiles, mobiles)
        slc2=np.ix_(np.asarray(range(nt)),immobiles, mobiles)
        massbalancecorrection=np.sum(dlnai_dlnwi[slc2]*wi0_immobiles[None,:,None],axis=1)
        #THFaktor[slc1]=dlnai_dlnwi[slc1]-massbalancecorrection[:,None,:]
        THFaktor=dlnai_dlnwi[slc1]-massbalancecorrection[:,None,:]
        THFaktor= interp1d(t,THFaktor,axis=0,bounds_error=False,fill_value=(THFaktor[0,:,:],THFaktor[-1,:,:]))
    
    
    #Set up the PDEs boundary,initial conditions and additional driving forces.
    #____________________________________
    #def default_mode():
    xinit=wvinit.flatten()
    dmuext=np.zeros((nTH,nz+1))
    wiB=wi8
    dwvdtB=np.zeros(nTH)
    def ode(t,x,THFaktor,dmuext,wiB,dwvdtB):
        """create generic ode fucnction for dwidt"""
        wv=np.ascontiguousarray(np.reshape(x,(nTH,nz+1))) 
        return drhodt(t,wv,THFaktor,mobiles,immobiles,Mi,D,allflux,swelling,rho,wi0,dmuext,wiB,dwvdtB)
    # Mechanical equation of state (MEOS)      
    if "EJ" in kwargs or "etaJ" in kwargs or "exponent" in kwargs: 
        from .relaxation import MEOS_mode
        xinit,ode=MEOS_mode(wvinit,ode,kwargs["EJ"],kwargs["etaJ"],kwargs["exponent"],Mi[mobiles],1/rho0i[mobiles])
    #_____________________________________
    if "witB" in kwargs: wiB=interp1d(t,kwargs['witB'],axis=0,bounds_error=False,fill_value=(kwargs['witB'][0],kwargs['witB'][-1]))

    def wrapode(t,x,ode,THFaktor,dmuext,wiB,dwvdtB):
        """evaluate time dependent functions and insert into the generic odes"""
        THFaktor=THFaktor(t)    if callable(THFaktor)   else THFaktor
        wiB=wiB(t)          if callable(wiB)      else wiB
        dwvdtB=dwvdtB(t,x)  if callable(dwvdtB)   else dwvdtB
        dxdt=ode(t,x,THFaktor,dmuext,wiB,dwvdtB)
        return dxdt.flatten()
    print("------------- Start diffusion modeling ----------------")
    start=time.time_ns()
    sol=solve_ivp(wrapode,(t[0],t[-1]),xinit,args=(ode,THFaktor,dmuext,wiB,dwvdtB),method="Radau",t_eval=t)
    end=time.time_ns()
    print("------------- Diffusion modeling took "+str((end-start)/1E9)+" seconds ----------------")
    if not sol["success"]: print("------------- Modeling failed the initial conditions are returned instead ----------------"); return wi0*np.ones((nc,nt)).T 
    x_sol=sol["y"] 
    nt=len(t)
    wt=np.zeros((nc,nt))
    wik=np.zeros((nt,nc,nz+1))
    wik=np.zeros((nt,nc,nz+1))
    for k in range(nt):
        wvk=np.reshape(x_sol[:(nz+1)*nTH,k],(nTH,nz+1))
        wik[k,:,:]=wiinit
        wik[k,mobiles,:]=wvk
        wik[k,:,-1]=wiB(t[k]) if callable(wiB) else wiB
        wik[k,immobiles,:]=(1-np.sum(wik[k,mobiles,:],axis=0))[None,:]*wi0_immobiles[:,None]
        wt[:,k]=np.sum(wik[k,:,:-1]/nz,axis=1)
    # Lt=1.wik/wt*L
    Lt=1
    # plt.plot(t,wik[:,1,-1])
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
        D[np.tril_indices_from(D,k=-1)]=Dvec
    return D

def Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=False,swelling=False,dlnai_dlnwi_fun=None,**kwargs):
    """iterates dlnai_dlnwi as a function of the concentration wi
    See also:
        diffusionpy.Diffusion_MS
    
    """
    nt=len(t)
    nc=len(wi0)
    if dlnai_dlnwi_fun is None : Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output,None,swelling,**kwargs)
    dlnai_dlnwi=np.stack([dlnai_dlnwi_fun(wi8*0.5+wi0*0.5)]*nt)
    wt_old=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,**kwargs)
    def wt_obj(wt_old):
        wt=wt_old.reshape((nt,nc))
        dlnai_dlnwi=np.asarray([dlnai_dlnwi_fun(wt[i,:]) for i in range(nt)])
        residual=(wt-Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,swelling=swelling,**kwargs)).flatten()/nt
        return residual
    wtopt=root(wt_obj,wt_old.flatten(),method="df-sane",options={"disp":True,"maxfev":30,"fatol":1E-6})["x"].reshape((nt,nc))

    if not full_output:
        return wtopt
    else: 
        dlnai_dlnwiopt=np.asarray([dlnai_dlnwi_fun(wtopt[i,:]) for i,val in enumerate(wtopt[:,0])])
        return Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=True,dlnai_dlnwi=dlnai_dlnwiopt,swelling=swelling,**kwargs)
    
def convert(x,M,axis=0):
    """convert fractions. e.g mass fractions into mole fractions"""
    y=np.empty(x.shape,M.dtype)
    np.copyto(y.T,M)
    return x*y/np.sum(x*y,axis=axis)
