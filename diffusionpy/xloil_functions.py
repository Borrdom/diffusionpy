try:
    import xloil as xlo
except ImportError:
    class xlo: 
        func = lambda f: f
        def Array(*args,**kwargs): return args[0]
    _has_xloil = False
else:
    _has_xloil = True
from .diffusion import Diffusion_MS,D_Matrix
import numpy as np
import time
from tkinter import Tk,filedialog,simpledialog
import xloil.pandas
import pandas as pd
from .surface import time_dep_surface
#from epcsaftpy import pcsaft,component,mixture
from .PCSAFT import lngi,vpure,dlnai_dlnxi
import fnmatch
# mix=component()+component()
# eos=pcsaft(mix)
from matplotlib import pyplot


@xlo.func
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
mobile:xlo.Array(bool,dims=1),swelling:bool=False,witB:xlo.Array(float,dims=2)=None,full_output:bool=False):   
    # if not _has_xloil: raise ImportError("xloil is required to do this.")
    kwargs={}
    if witB is not None: kwargs.update({"witB":witB.copy()})
    if full_output: 
        _,wtz,zvec,Lt=Diffusion_MS(t.copy(),L,Dvec.copy(),w0.copy(),w8.copy(),mobile,True,None,swelling,**kwargs)
        nz=len(zvec)
        zt=Lt[None,:]*zvec[:,None]/L
        res=zt.flatten()[::-1,None]
        for i in range(len(w0)):
            res=np.hstack((res,wtz[:,i,:].flatten()[:,None]))
        return res
    else: 
        return Diffusion_MS(t.copy(),L,Dvec.copy(),w0.copy(),w8.copy(),mobile,**kwargs)[0]
@xlo.func
def Diffusion_MS_iter_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
mobile:xlo.Array(bool,dims=1),swelling:bool=False,witB:xlo.Array(float,dims=2)=None,T:float=298.15,p:float=1E5,
pure:xlo.Array(object,dims=2)=np.asarray([[]]),kij:xlo.Array(object,dims=2)=np.asarray([[]]),maxit:int=10,full_output:bool=False):
    Mw=pure[1:,3].astype(float)
    mi=pure[1:,4].astype(float)
    sigi=pure[1:,6].astype(float)
    ui=pure[1:,7].astype(float)
    eAiBi=pure[1:,8].astype(float)
    kAiBi=pure[1:,9].astype(float)
    Na=pure[1:,10].astype(float)
    nc=len(Mw)
    kij1=np.char.replace(kij[1:,2].astype(str),",",".").astype(float)
    kijhb1=np.char.replace(kij[1:,3].astype(str),",",".").astype(float)
    # kij=D_Matrix(kij1,nc)
    # kijhb=D_Matrix(kij1hb,nc)
    vpures=vpure(p,T,mi,sigi,ui,eAiBi,kAiBi,Na)
    par={"mi":mi,
    "Mi":Mw,
    "si":sigi,
    "ui":ui,
    "eAi":eAiBi,
    "kAi":kAiBi,
    "NAi":Na,
    "kij":kij1,
    "kijA":kijhb1,
    "vpure":vpures}
    dlnai_dlnwi_fun=lambda wi: dlnai_dlnxi(T,wi,**par)
    kwargs={}
    if witB is not None: kwargs.update({"witB":witB.copy()})
    kwargs.update({"method":"fixedpoint"})
    kwargs.update({"maxiter":7})
    if full_output: 
        _,wtz,zvec,Lt=Diffusion_MS(t.copy(),L,Dvec.copy(),w0.copy(),w8.copy(),mobile,saftpar=par,**kwargs)
        nz=len(zvec)
        zt=Lt[None,:]*zvec[:,None]/L
        res=zt.flatten()[::-1,None]
        for i in range(len(w0)):
            res=np.hstack((res,wtz[:,i,:].flatten()[:,None]))
        return res
    else: 
        return Diffusion_MS(t.copy(),L,Dvec.copy(),w0.copy(),w8.copy(),mobile,saftpar=par,**kwargs)[0]


@xlo.func
def time_dep_surface_xloil(t:xlo.Array(float,dims=1),wi0:xlo.Array(float,dims=1),wi8:xlo.Array(float,dims=1),mobile:xlo.Array(bool,dims=1),taui:xlo.Array(float,dims=1),lngi_t:xlo.Array(float,dims=2)=None):
    return time_dep_surface(t.copy(),wi0.copy(),wi8.copy(),mobile.copy(),taui.copy(),lngi_t)

@xlo.func
def gradient(x:xlo.Array(float,dims=1),y:xlo.Array(float,dims=1)):
    return np.gradient(x,y)


@xlo.func
def reduce_points(x,n:int):
    return x[::n]

@xlo.func
def interp1d(x,xp,fp):
    return np.interp(x.flatten(),xp.flatten(),fp.flatten())
