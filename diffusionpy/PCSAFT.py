import numpy as np
from scipy.optimize import root
from  .ares import ares 
#@njit(cache=True)
def eta_iter(p,T,xi,**kwargs):
    """solve the density mich yiels a given pressure p"""
    xi=np.ascontiguousarray(xi)
    def Z_obj(p,T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA):
        kB = 1.380649e-23
        di=si*(1.-0.12*np.exp(-3*ui/T))
        rho=6/np.pi*eta*(np.sum(mi*xi*di**3))**-1
        rhobar=rho*(10.**10)**3
        Zp=p/(rhobar*kB*T)
        _,_,Z1=ares(T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA)
        return (Zp-Z1.real)
    def eta_roots(fun,p,T,eta,xi,**kwargs):
        f=fun(p,T,eta,xi,**kwargs)#.reshape(n)
        tol=1E-8
        it=50
        h = tol
        deta = h
        J= (fun(p,T,eta+deta,xi,**kwargs)-f)/h
        for i in range(it):
            if np.abs(f)<tol:
                return eta 
            s=-1.*f/J
            eta+=s
            df=fun(p,T,eta,xi,**kwargs)-f
            J+=(df-J*s)/s
            f+=df
        return eta
    eta0=0.45
    return eta_roots(Z_obj,p,T,eta0,xi,**kwargs)

def vpure(p,T,mi,si,ui,eAi,kAi,NAi,**kwargs):
    """solve the density mich yiels a given pressure p"""
    etapures=[]
    for i in range(len(mi)):
        purepar={'mi':np.asarray([mi[i]]),
        'si':np.asarray([si[i]]),
        'ui':np.asarray([ui[i]]),
        'eAi':np.asarray([eAi[i]]),
        'kAi':np.asarray([kAi[i]]),
        'NAi':np.asarray([NAi[i]]),
        'kij':np.zeros(10),
        'kijA':np.zeros(10)}
        x=eta_iter(p,T,np.asarray([1.]),**purepar)
        etapures.append(x)
    etapures=np.asarray(etapures)
    di=si*(1.-0.12*np.exp(-3*ui/T))
    NA = 6.0221407e23
    vmol=np.pi/6/etapures*mi*di**3/(10.**10)**3*NA
    return vmol

#@njit(cache=True)
# @njit(['Tuple((f8[::1], f8[::1], f8[::1],f8[::1]))(f8,f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])',
# 'Tuple((c16[::1], c16[::1], c16[::1],c16[::1]))(f8,c16[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])'],cache=True)
def SAFTSAC(T,xi,mi,si,ui,eAi,kAi,NAi,vpure,Mi,kij,kijA):
    """Calculate the log of the activity coefficients via the SAFT-SAC approximation
    Args:
        T (float): temperature
        xi (array_like): mole/mass fraction. Becomes the mass fraction when the molar mass Mi is not None
        mi (array_like): segment number
        si (array_like): segment diameter
        ui (array_like): dispersion energy
        eAi (array_like): association energy
        kAi (array_like): association volume
        NAi (array_like): association sites (only symmetric)
        vpure (array_like): pure component molar volumes
        Mi (array_like, optional): Molar mass. Calculates properties on a mass basis when given. Defaults to None.
        kij (array_like, optional): Matrix of binary interaction parameters for dispersion . Defaults to np.asarray([[0.]]).
        kijA (array_like, optional): Matrix of binary interaction parameters for association Defaults to np.asarray([[0.]]).
    Returns:
        array_like: vector of activity coefficients
    """
    NA = 6.0221407e23
    xi=np.ascontiguousarray(xi)
    #vpfracNET=(1-ksw*RH**2)/xi[0]
    #vmol=v0pNE/vpfracNET
    wi=xi
    #if Mi is not None: xi=xi/Mi/np.sum(xi/Mi)
    xi=xi/Mi/np.sum(xi/Mi)
    vmol=np.sum(vpure*xi) #if "vmol" not in kwargs else kwargs["vmol"]
    vpfrac=vpure/vmol
    # vpfrac=mi/np.sum(mi*xi)
    di=si*(1.-0.12*np.exp(-3*ui/T))
    eta=np.pi/6*np.sum(mi*xi.real*di**3)/vmol/(10.**10)**3*NA
    etapure=np.pi/6*mi*di**3/vpure/(10.**10)**3*NA
    lngi_id=np.log(vpfrac)+1.-vpfrac
    arespures=np.asarray([ares(T,val,np.asarray([1.]),np.asarray([mi[i]]),np.asarray([si[i]]),np.asarray([ui[i]]),np.asarray([eAi[i]]),np.asarray([kAi[i]]),np.asarray([NAi[i]]),np.zeros(10),np.zeros(10))[0] for i,val in enumerate(etapure)])
    _,mures,Z1=ares(T,eta,xi,mi,si,ui,eAi,kAi,NAi,kij,kijA)
    lngi_res=mures-arespures
    lngi_p=-vpfrac*(Z1-1) #if "NETGP" not in kwargs else 0
    #with np.errstate(divide='ignore',invalid='ignore'):
    with np.errstate(divide='ignore',invalid='print'):
        lngi_wx=np.nan_to_num(np.log(xi/wi),0)
    return lngi_id,lngi_res,lngi_p,lngi_wx

# @njit(['f8[::1](f8,f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])',
# 'c16[::1](f8,c16[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])'],cache=True)
def lngi(T,xi,mi,si,ui,eAi,kAi,NAi,vpure,Mi,kij,kijA,**kwargs):
    lngi_id,lngi_res,lngi_p,lngi_wx=SAFTSAC(T,xi,mi,si,ui,eAi,kAi,NAi,vpure,Mi,kij,kijA)
    return lngi_id+lngi_res+lngi_p+lngi_wx

# @njit('f8[:,::1](f8,f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])',cache=True)
def dlnai_dlnxi(T,xi,**kwargs):
    """Generate the derivatives of the mole fraction with concentration

    Args:
        T (float): temperature
        xi (array_like): mole/mass fraction. Becomes the mass fraction when the molar mass Mi is not None
        par (dic) : dictionary containg pc-saft parameters
    Return:
        array_like: martrix of derivatives of the mole fraction with concentration
    """
    nc = len(xi)
    h = 1E-26
    df = np.zeros((nc, nc))
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'c16')
        dx[i] = h * 1j
        out =  lngi(T,xi+dx,**kwargs)
        df[i] = out.imag/h
    return df.T*xi+np.eye(nc)

# @njit('f8[:,:,:,::1](f8,f8[:,:,::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[::1],f8[:,:],f8[:,:])',cache=True)
def dlnai_dlnxi_loop(T,xi,**kwargs):
    nt,nz,nc = xi.shape
    dlnai_dlnxi_vec=np.zeros((nt,nz,nc,nc))
    for i in range(nt):
        for j in range(nz):
            dlnai_dlnxi_vec[i,j,:,:]=dlnai_dlnxi(T,np.ascontiguousarray(xi[i,j,:]),**kwargs)
    return dlnai_dlnxi_vec

def NETVLE(T,wi,v0p,mobile,polymer,ksw,mi,sigi,ui,epsAiBi,kapi,N,vpures,Mi,kij,kijA,n=2):
    # vp=np.zeros(np.sum(polymers))
    
    mobiles=np.where(mobile)[0] #if not allflux else np.arange(0,nc-1,dtype=np.int64)
    immobiles=np.where(~mobile)[0] #if not allflux else np.asarray([nc-1],dtype=np.int64)
    widry=np.zeros_like(wi)
    widry[immobiles]=wi[immobiles]/np.sum(wi[immobiles])
    nc=len(Mi)
    kswmat=np.zeros((nc,nc))
    polymers=np.where(polymer)[0]
    v0NET=vpures*1000./Mi
    v0NET[polymers]=v0p
    Mdry=(np.sum(widry/Mi))**-1
    xidry=widry/Mi*Mdry
    v0dry=np.sum(v0NET*widry)
    kswmat[mobiles,polymers]=ksw
    vidry=v0NET/v0dry*widry
    ksws=(kswmat@vidry)[mobiles]
    v0moldry=v0dry*Mdry/1000.
    def res(RS):
        # RS=np.fmin(np.fmax(RS,1E-8),1)
        xi=wi/Mi/np.sum(wi/Mi,axis=0)
        xw=xi[mobiles]
        ww=wi[mobiles]
        vmol=v0moldry/(1-np.sum(ksws*RS**n))*(1-np.sum(xw))
        vmoltrick=(vmol-np.sum(xw*vpures[mobiles]))/(1-np.sum(xw))
        # v=v0dry/(1-np.sum(ksws*RS**2))*(1-np.sum(ww))
        # vtr=(v-np.sum(ww*v0NET[mobiles]))/(1-np.sum(ww))*widry[immobiles]
        vpures2=vpures.copy()
        vpures2[immobiles]=np.fmax(vmoltrick,1E-12)
        # vpures2[immobiles]=vtr*Mi[immobiles]/1000.
        # vpures2[immobiles]=(vmol-np.sum(xw*vpures[mobiles]))/(1-np.sum(xw))*xidry[immobiles]
        lngid,lngres,_,lngw=SAFTSAC(T,wi,mi,sigi,ui,epsAiBi,kapi,N,vpures2,Mi,kij,kijA)
        logRS=lngid[mobiles]+lngres[mobiles]+lngw[mobiles]+np.log(wi[mobiles])
        return logRS-np.log(RS)
    re=root(res,wi[mobiles]/2,method='hybr')
    RS=re["x"]
    return RS

def supersaturation(T,xi,mi,si,ui,eAi,kAi,NAi,vpure,Mi,kij,kijA,deltaHSL,TSL,cpSL):
    R=8.3145
    lnaiSLE=-deltaHSL/(R*T)*(1-T/TSL)+cpSL/R*(TSL/T-1-np.log(TSL/T))
    lnai=lngi(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij, kijA) + np.log(xi)
    return lnai-lnaiSLE


# class mixture:
#     def __init__(self,p,T,mi,si,ui,eAi,kAi,NAi,Mi,kij,kijA):
#         self.mi=mi
#         self.si=si
#         self.ui=ui
#         self.eAi=eAi
#         self.kAi=kAi
#         self.NAi=NAi
#         self.Mi=Mi
#         self.kij=kij
#         self.kijA=kijA
#         self.vpure=vpure(p,T, **vars(self))