
import numpy as np
from scipy.interpolate import interp1d
from diffusionpy import D_Matrix,Diffusion_MS_iter,Diffusion_MS,time_dep_surface_cryst,lngi,vpure,dlnai_dlnxi,origin_like
import matplotlib.pyplot as plt

crystallize=np.asarray([False,False,True])
mobile=np.asarray([True,False,False])
mobiles=np.where(mobile)[0]
immobiles=np.where(~mobile)[0]
deltaHSL=np.asarray([31500.])
TSL=np.asarray([429.47])
cpSL=np.asarray([87.44])


texp=np.asarray([15,30,45,60,90,120])

## Dl10
relapi10=np.asarray([10.50583658,19.84435798,28.40466926,38.13229572,52.14007782,62.25680934])
relpoly10=np.asarray([8.62745098,19.21568627,28.62745098,36.8627451,49.01960784,60.78431373])

## Dl20
relapi20=np.asarray([5.836575875,6.614785992,13.61867704,19.06614786,27.62645914,36.96498054])
relpoly20=np.asarray([5.098039216,8.235294118,16.07843137,22.74509804,32.94117647,42.74509804])

## Dl30
relapi30=np.asarray([1.945525292,1.945525292,3.112840467,3.891050584,5.447470817,7.003891051])
relpoly30=np.asarray([1.176470588,1.568627451,3.529411765,5.882352941,6.274509804,9.019607843]) 

DAPI=np.asarray([6.6E-17])
sigma=np.asarray([1E-02])
kt=np.asarray([5.1E-12])
g=np.asarray([3.2])
Mi=np.asarray([18.015,65000.,230.26])
rho0i=np.asarray([997.,1180.,1320.])

nc=3
wv0=0.0001
dl0=0.1
wv8=0.99
def limits(dl0,wv0,wv8):
       wi0=np.asarray([wv0,(1-wv0)*(1-dl0),(1-wv0)*dl0])
       wi8=np.asarray([wv8,(1-wv8)*(1-dl0),(1-wv8)*dl0])
       return wi0,wi8
T=298.15
p=1E5


kij=D_Matrix(np.asarray([-0.128,0.00648,-0.0574]),nc)
par={"mi":np.asarray([1.20469,2420.99, 8.105152]),
"si": np.asarray([2.797059952,2.947, 2.939]),
"ui" :np.asarray([353.95,205.27, 229.45]),
"eAi" :np.asarray([2425.67,0., 934.2]),
"kAi":np.asarray([0.04509,0.02, 0.02]),
"NAi":np.asarray([1.,653., 2.]),
"Mi": Mi,
"kij":kij}

vpures=vpure(p,T,**par)

par["vpure"]=vpures
lngi_fun=lambda wi: lngi(T,wi,**par)
R=8.3145
lnaiSLE=-deltaHSL/(R*T)*(1-T/TSL)+cpSL/R*(TSL/T-1-np.log(TSL/T))
nt=30
t=np.linspace(0,texp[-1],nt)*60
Dvec=np.asarray([1E-14,1E-15,1E-15])
L=3E-5

def iterate(wi0,wi8):
       wt,wtz,zvec,Lt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,full_output=True)
       for i in range(5):
              lngi_tz=np.asarray([[lngi_fun(col) for col in row.T] for row in wtz])
              wt,wtz,zvec,Lt,_,_=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,deltaHSL=deltaHSL,TSL=TSL,cpSL=cpSL,crystallize=crystallize,sigma=sigma,DAPI=DAPI,kt=kt,g=g,full_output=True,lngi_tz=lngi_tz)
              notreleased=wt/wi0
              release=(1-notreleased)
       return release

dl010=0.1
dl020=0.2
dl030=0.3

wi010,wi810=limits(dl010,wv0,wv8)
wi020,wi820=limits(dl020,wv0,wv8)
wi030,wi830=limits(dl030,wv0,wv8)

release10=iterate(wi010,wi810)
release20=iterate(wi020,wi820)
release30=iterate(wi030,wi830)
# plt.plot(t/60,witB.T[:,0]"bo")

# plt.plot(t/60,release[:,1]*100,label="Non-ideal and wiB(t)")
# plt.plot(t/60,alpha*100,'r',label="Non-ideal and wiB(t)")
fig,ax=origin_like.subplots()
ax.set_axisbelow(True)
origin_like.plot(ax,t/60,release10[:,2]*100,"g-")
origin_like.plot(ax,t/60,release20[:,2]*100,"c-")
origin_like.plot(ax,t/60,release30[:,2]*100,"-r")

origin_like.plot(ax,texp,relapi10,"go")
origin_like.plot(ax,texp,relapi20,"cs")
origin_like.plot(ax,texp,relapi30,"^r")

ax.set_xlabel("t/min")
ax.set_ylabel("API release/%")
origin_like.set_ticks(ax,0,120,0,100)
plt.ylim(0,100)
plt.xlim(0,130)


plt.show()
# import pandas as pd
# pd.DataFrame((tmin,ww)).T.to_clipboard(excel=True, sep=None, index=False, header=None)


