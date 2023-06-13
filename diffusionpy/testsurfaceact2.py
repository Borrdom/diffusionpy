import matplotlib.pyplot as plt
import numpy as np
from .Stefan_Maxwell_segmental import Diffusion_MS_iter,D_Matrix,Diffusion_MS
from .PyCSAFT_nue import vpure
from .surface_activity import time_dep_surface2
import matplotlib.pyplot as plt

nt=500
t=np.linspace(0,429.8518201,nt)
Dvec=np.asarray([1E-6,2.3E-11,1.7E-11])
nc=3
L=2E-5
wi0=np.asarray([0.333333333,0.333333333,0.333333333])
wi8=np.asarray([0.00001,0.127606346,0.872393654])
Mi=np.asarray([32.042,92.142,90000.])
mobile=np.asarray([True,True,False])


from .PyCSAFT_nue import dlnai_dlnxi,vpure,lngi
T=298.15
p=1E5
Dvec=np.asarray([1E-6,2.3E-10,1.7E-10])
kij=D_Matrix(np.asarray([0.029,-0.05855362,0.027776682]),nc)
par={"mi" :np.asarray([1.5255, 2.8149, 2889.9]),
"ui" : np.asarray([188.9, 285.69, 204.65]),
"si" : np.asarray([3.2300, 3.7169, 3.3972]),
"kAi": np.asarray([0.035176, 0., 0.]),
"eAi": np.asarray([2899.5, 0., 0.]),
"NAi": np.asarray([1., 0., 1047.]),
"Mi" : np.asarray([32.042,  92.142, 90000.]),
"kij": kij}
vpures=vpure(p,T,**par)
par["vpure"]=vpures
dlnai_dlnwi=np.stack([dlnai_dlnxi(T,wi8*0.5+wi0*0.5,**par)]*nt)
wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par)
plt.plot(t,wt[:,0])
plt.plot(t,wt[:,1])
plt.plot(t,wt[:,2])

taui=np.asarray([1,70])
from .surface_activity import time_dep_surface2
witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=None)
wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par,witB=witB)
# wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]))
# wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]),par=par)
plt.plot(t,wt[:,0])
plt.plot(t,wt[:,1])
plt.plot(t,wt[:,2])

witB=time_dep_surface2(t,wi0,wi8,mobile,taui,par=par)
wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,T=T,par=par,witB=witB)
# wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]))
# wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,dlnai_dlnwi=dlnai_dlnwi,taui=np.asarray([1,70]),par=par)
plt.plot(t,wt[:,0])
plt.plot(t,wt[:,1])
plt.plot(t,wt[:,2])
plt.savefig('filename.png', format='png',  transparent=True)
plt.show()