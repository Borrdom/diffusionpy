# %%
import numpy as np
from scipy.interpolate import interp1d
from diffusionpy import D_Matrix,Diffusion_MS_iter,Diffusion_MS,time_dep_surface_cryst,lngi,vpure,dlnai_dlnxi
import matplotlib.pyplot as plt

# %%
crystallize=np.asarray([False,False,True])
mobile=np.asarray([True,False,False])
mobiles=np.where(mobile)[0]
immobiles=np.where(~mobile)[0]
deltaHSL=np.asarray([31500.])
TSL=np.asarray([429.47])
cpSL=np.asarray([87.44])

# %% [markdown]
# The experimental data was taken from https://doi.org/10.3390/pharmaceutics15051539
# 
#  # Dl10 37C
# ![image-2.png](attachment:image-2.png) 0 min
# ![image-3.png](attachment:image-3.png) 10 min
# ![image-7.png](attachment:image-7.png) 20 min
# ![image-4.png](attachment:image-4.png) 30 min
# ![image-5.png](attachment:image-5.png) 40 min
# ![image.png](attachment:image.png)     50 min
# ![image-6.png](attachment:image-6.png) 60 min

# %% [markdown]
# The experimental data was taken from https://doi.org/10.3390/pharmaceutics15051539
# 
#  # Dl30 37C
#  ![image.png](attachment:image.png)      0 min
#  ![image-2.png](attachment:image-2.png) 10 min
#  ![image-3.png](attachment:image-3.png) 20 min
#  ![image-4.png](attachment:image-4.png) 30 min
#  ![image-5.png](attachment:image-5.png) 40 min
#  ![image-6.png](attachment:image-6.png) 50 min
#  ![image-7.png](attachment:image-7.png) 60 min

# %%
texp=np.asarray([15,
30,
45,
60,
90,
120])

## Dl10
relapi10=np.asarray([10.50583658,
19.84435798,
28.40466926,
38.13229572,
52.14007782,
62.25680934])

relpoly10=np.asarray([8.62745098,
19.21568627,
28.62745098,
36.8627451,
49.01960784,
60.78431373])

## Dl20

relapi20=np.asarray([5.836575875,
6.614785992,
13.61867704,
19.06614786,
27.62645914,
36.96498054])

relpoly20=np.asarray([5.098039216,
8.235294118,
16.07843137,
22.74509804,
32.94117647,
42.74509804])


## Dl30

relapi30=np.asarray([1.945525292,
1.945525292,
3.112840467,
3.891050584,
5.447470817,
7.003891051])

relpoly30=np.asarray([1.176470588,
1.568627451,
3.529411765,
5.882352941,
6.274509804,
9.019607843])

# %% [markdown]
# 

# %%
data=np.asarray([[0.,	0.27288],
       [0.05   , 0.23445],
       [0.1    , 0.19755],
       [0.15   , 0.17761],
       [0.2    , 0.16723],
       [0.25   , 0.15672],
       [0.3    , 0.14609],
       [0.35   , 0.13535],
       [0.45   , 0.11365],
       [0.5    , 0.10276],
       [0.55   , 0.09191],
       [0.6    , 0.0812 ],
       [0.65   , 0.07075],
       [0.7    , 0.0607 ],
       [0.75   , 0.05123],
       [0.8    , 0.04251],
       [0.85   , 0.03468],
       [0.9    , 0.03224],
       [0.95   , 0.0312 ],
       [1.     , 0.03033]])
ww=data[:,1] #np.asarray([0.27087,0.22302, 0.13792, 0.09208, 0.06118])
dl=data[:,0] #np.asarray([0,0.1 , 0.3 , 0.5 , 0.68])
DAPI=np.asarray([6.6E-17])
sigma=np.asarray([1E-02])
kt=np.asarray([5.1E-12])
g=np.asarray([3.2])
Mi=np.asarray([18.015,65000.,230.26])
rho0i=np.asarray([997.,1180.,1320.])

# %%
wv_fun=interp1d(dl,ww,bounds_error=False,fill_value=(np.max(ww),np.min(ww)))

# %%
nc=3
wv0=0.0001
dl0=0.3
wi0=np.asarray([wv0,(1-wv0)*(1-dl0),(1-wv0)*dl0])
wv8=0.99#wv_fun(dl0)
wi8=np.asarray([wv8,(1-wv8)*(1-dl0),(1-wv8)*dl0])
T=298.15
p=1E5

# %%
kij=D_Matrix(np.asarray([-0.128,0.00648,-0.0574]),nc)
par={"mi":np.asarray([1.20469,2420.99, 8.105152]),
"si": np.asarray([2.797059952,2.947, 2.939]),
"ui" :np.asarray([353.95,205.27, 229.45]),
"eAi" :np.asarray([2425.67,0., 934.2]),
"kAi":np.asarray([0.04509,0.02, 0.02]),
"NAi":np.asarray([1.,653., 2.]),
"Mi": Mi,
"kij":kij}

# %%
vpures=vpure(p,T,**par)

# %%
par["vpure"]=vpures
lngi_fun=lambda wi: lngi(T,wi,**par)
R=8.3145
lnaiSLE=-deltaHSL/(R*T)*(1-T/TSL)+cpSL/R*(TSL/T-1-np.log(TSL/T))
lnai0=lngi_fun(wi0)+np.log(wi0)

# lngi_fun=lambda wi: 
# %% [markdown] 
# # Ives–Pienvichitr Power Law Relation https://doi.org/10.1007/s11242-018-1086-2
# # Alternativ Krishna https://doi.org/10.1016/S0009-2509(96)00458-7
# #und dann Knudsen Diffusion. Die Porengröße und ide porosität ändert sich mit alpha
# 
# eta=1.5# hängt von der Toriszität ab
# 
# #http://www.idaea.csic.es/mhetscale/wp-content/uploads/2019/06/DiffusionReviewPreprint.pdf
# # Deff=amorphizität/tau*D wobei tau=1/amorphizität**0.5
# Lambda=amorphizität**eta  # referenzamorphizitä 1 dadurch voller diffusionkoeffizient effektiver diffusionskoeffizient

# %%
nt=30
t=np.linspace(0,texp[-1],nt)*60
witB=np.ones((nc,nt))*wi0[:,None]
lngi_t=interp1d(t,np.asarray([lngi_fun(wi0)]*nt),axis=0)
for i in range(5):
       witB,alpha,r=time_dep_surface_cryst(t,np.asarray([True,False,False]),np.asarray([False,True,True]),crystallize,wi0,wi8,rho0i,Mi,DAPI,sigma,kt,g,deltaHSL,TSL,cpSL,lngi_t,wv_fun)
       lngi_t=interp1d(t,np.asarray([lngi_fun(wi) for wi in witB.T]),axis=0)
       # plt.plot(t,alpha)
tmin=t/60
# plt.show()
Dvec=np.asarray([1E-14,1E-15,1E-15])
L=3E-5



# dlnai_dlnwi_fun=lambda wi: dlnai_dlnxi(T,wi,**par)
# dlnai_dlnwi=np.stack([dlnai_dlnwi_fun(wi8*0.5+wi0*0.5)]*nt)
# dlnai_dlnwi=np.stack([np.eye(3,3)]*nt)
# # amorphizität=(1-alpha)[:,None,None]
# Lambda=1-np.exp(-amorphizität/0.3)

# eta=1.5
# Lambda=amorphizität**eta  # referenzamorphizitä 1 dadurch voller diffusionkoeffizient effektiver diffusionskoeffizient

# wt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,dlnai_dlnwi=dlnai_dlnwi*Lambda)
wt,wtz,zvec,Lt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=True)

for i in range(2):
       lngi_tz=np.asarray([[lngi_fun(col) for col in row.T] for row in wtz])
       wt,wtz,zvec,Lt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,deltaHSL=deltaHSL,TSL=TSL,cpSL=cpSL,crystallize=crystallize,sigma=sigma,DAPI=DAPI,kt=kt,g=g,full_output=True,lngi_tz=lngi_tz)
# lngi_tz=np.asarray([[lngi_fun(col) for col in row.T] for row in wtz])
# wt,wtz,zvec,Lt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,deltaHSL=deltaHSL,TSL=TSL,cpSL=cpSL,crystallize=crystallize,sigma=sigma,DAPI=DAPI,kt=kt,g=g,full_output=True,lngi_tz=lngi_tz)
# wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,witB=witB.T,dlnai_dlnwi_fun=dlnai_dlnwi_fun)
       notreleased=wt/wi0
       release=(1-notreleased)

       wtid=wt
       # plt.plot(t/60,witB.T[:,0],"bo")
       plt.plot(t/60,release[:,2]*100,label="Non-ideal and wiB(t)")
       # plt.plot(t/60,release[:,1]*100,label="Non-ideal and wiB(t)")
# plt.plot(t/60,alpha*100,'r',label="Non-ideal and wiB(t)")

plt.plot(texp,relapi10,'bx',label="Non-ideal and wiB(t)")
plt.plot(texp,relapi20,'bx',label="Non-ideal and wiB(t)")
plt.plot(texp,relapi30,'bx',label="Non-ideal and wiB(t)")

plt.xlim(0,130)
plt.ylim(0,100)
plt.legend("")
plt.show()
import pandas as pd
pd.DataFrame((tmin,ww)).T.to_clipboard(excel=True, sep=None, index=False, header=None)


