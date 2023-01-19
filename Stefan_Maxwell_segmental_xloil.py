import xloil as xlo
from Stefan_Maxwell_segmental import Diffusion_MS,Diffusion1D

@xlo.func()
def Diffusion_MS_xloil(t:xlo.Array(float,dims=1),L:float,Dvec:xlo.Array(float,dims=1),w0:xlo.Array(float,dims=1),w8:xlo.Array(float,dims=1),Mi:xlo.Array(float,dims=1),volatile:xlo.Array(bool,dims=1),full_output=False):
    return Diffusion_MS(t,L,Dvec,w0,w8,Mi,volatile,full_output)

@xlo.func()
def Diffusion1D_xloil(t:xlo.Array(float,dims=1),L0:float,Ds:float,ws0:float,ws8:float):
    return Diffusion1D(t,L0,Ds,ws0,ws8)

