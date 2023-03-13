import xloil as xlo
from .Stefan_Maxwell_segmental import Diffusion_MS,Diffusion1D
from .crank_and_other import crank,BHX
from .DasDennisSpacing import DasDennis
import numpy as np
import casadi as cs
import time

@xlo.func
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
volatile:xlo.Array(bool,dims=1),full_output:bool,Gammai:xlo.Array(float,dims=2)=None,swelling:bool=False,taui:xlo.Array(float,dims=1)=None):   
    return Diffusion_MS(t,L,Dvec,w0,w8,Mi,volatile,full_output,Gammai,swelling,taui)

@xlo.func
def Diffusion1D_xloil(t:xlo.Array(float,dims=1),L0:float,Ds:float,ws0:float,ws8:float):
    return Diffusion1D(t,L0,Ds,ws0,ws8)
    
@xlo.func
def gradient(x:xlo.Array(float,dims=1),y:xlo.Array(float,dims=1)):
    return np.gradient(x,y)

@xlo.func
def reduce_points(x,n:int):
    return x[::n]
@xlo.func
def crank_xl(t,L0,Ds,ws0,ws8):
    Xs0=ws0/(1-ws0)
    Xs8=ws8/(1-ws8)
    mt08=crank(t,Ds/L0**2,1)
    Xs=mt08*(Xs8-Xs0)+Xs0
    ws=Xs/(1+Xs)
    return ws
@xlo.func
def interp1d(x,xp,fp):
    return np.interp(x.flatten(),xp.flatten(),fp.flatten())

@xlo.func
def BHX_xloil(t:xlo.Array(float,dims=1),kf:float,kr:float,ws0:float,ws8:float,mfinfty:float,mrinfty:float):
    Xs0=ws0/(1-ws0)
    Xs8=ws8/(1-ws8)
    Xs=BHX(t,kf,np.asarray([kr]),1,Xs0,Xs8,mfinfty,np.asarray([mrinfty]))
    ws=Xs/(1+Xs)
    return ws

@xlo.func
def DasDennis_xloil(p:int,dim:int):
    return DasDennis(p,dim).T

@xlo.func
def add_custom(a,b):
    return a+b