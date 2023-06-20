from __future__ import division, print_function, absolute_import
import numpy as np
from .Stefan_Maxwell_segmental import Diffusion_MS,D_Matrix,Diffusion_MS_iter,convert
from .surface_activity import time_dep_surface
from .xloil_functions import Diffusion_MS_xloil,reduce_points,crank_xl,interp1d,BHX_xloil,DasDennis_xloil,Diffusion_MS_iter_xloil
from .PyCSAFT_nue import lngi,eta_iter,ares,lnphi_TP,vpure,dlnai_dlnxi


T=298.15
p=1E5
npoint=2
#Water
#XAi0 =0.16 und 0.6

par={"mi":np.asarray([1.20469,2.38269789]),
 "si":np.asarray([2.797059952,3.1771]),
 "ui":np.asarray([353.95,198.24]),
 "eAi":np.asarray([2425.67,2653.4]),
 "kAi":np.asarray([0.04509,0.032384]),
 "NAi":np.asarray([1.,1.])}
x1=np.linspace(0,1,npoint)
x2=1-x1
xi=np.vstack((x1,x2))
vpures=vpure(p,T,**par)
par["vpure"]=vpures
lngiammai=np.asarray([lngi(T,xi[:,i],**par).flatten() for i,val in enumerate(xi[0,:])])
Gammai=np.asarray([dlnai_dlnxi(T,xi[:,i],**par).flatten() for i,val in enumerate(xi[0,:])])