import xloil as xlo
from PyFusion.Stefan_Maxwell_segmental import Diffusion_MS,Diffusion1D
import numpy as np
import casadi as cs
import time

@xlo.func
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),
volatile:xlo.Array(bool,dims=1),full_output:bool,Gammai:xlo.Array(float,dims=2)):   
    return Diffusion_MS(t,L,Dvec,w0,w8,Mi,volatile,full_output,Gammai)

@xlo.func
def Diffusion1D_xloil(t:xlo.Array(float,dims=1),L0:float,Ds:float,ws0:float,ws8:float):
    return Diffusion1D(t,L0,Ds,ws0,ws8)
    
@xlo.func
def gradient(x:xlo.Array(float,dims=1),y:xlo.Array(float,dims=1)):
    return np.gradient(x,y)