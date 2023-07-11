# %% [markdown]
# # Stefan-Maxwell diffusion involving three components
# 

# %% [markdown]
# This example demonstrates the modeling of the multicomponent Stefan-Maxwell model involving three componentt.
# First we import the Stefan-Maxwell diffusion module

# %%
import numpy as np
from numpy import array
from diffusionpy import Diffusion_MS,D_Matrix,Diffusion_MS_iter,vpure,dlnai_dlnxi,lngi,circular
import matplotlib.pyplot as plt


# %% [markdown]
# experimental data

# %%

swelling=array([[0.       ,    0.    ],
       [  9.24861,   0.895  ],
       [ 11.18739,   1.08754],
       [ 16.86526,   1.2824 ],
       [ 32.37552,   1.54222],
       [ 63.25758,   2.15696],
       [121.14411,   2.80418],
       [152.02616,   3.06168],
       [182.76973,   3.28902],
       [211.713  ,   3.51403],
       [242.59505,   3.8388 ],
       [271.53832,   4.00119],
       [302.42037,   4.25869],
       [333.30242,   4.38859],
       [360.3069 ,   4.61593],
       [420.13222,   5.06829],
       [451.01427,   5.13324],
       [483.83511,   5.35826]])
texp=swelling[:,0]
thickness=swelling[:,1]
msol=500000
mges0=125
ww0=0.01
dl0=0.4
dl8=0.4
wi0=np.asarray([(1-dl0)*(1-ww0),dl0*(1-ww0),ww0])
# release=wL1D03*msol/m0
# notreleased=1-release
# wexp=notreleased*wi0
# wexp[:,2]=1-wexp[:,0]-wexp[:,1]

wexp=thickness



# %% [markdown]
# We want to describe the diffusion of water into an ASD film. The ASD-water mixture is a ternary system. First we define starting and equilibrium weight fractions.

# %%
nc=3
L=6/1000
ww8=0.95
wi8=np.asarray([(1-dl8)*(1-ww8),dl8*(1-ww8),ww8])
Mi=np.asarray([86000,554.5,18.015])
T=298.15
p=1E5

kij=D_Matrix(np.asarray([-0.0621,-0.156,-0.025]),nc)
par={"mi":np.asarray([2420.99, 14.283,1.20469]),
"si": np.asarray([2.947, 3.535,2.797059952]),
"ui" :np.asarray([205.27, 262.79,353.95]),
"eAi" :np.asarray([0., 886.4,2425.67]),
"kAi":np.asarray([0.02, 0.02,0.04509]),
"NAi":np.asarray([653., 3.,1.]),
"Mi": Mi,
"kij":kij}

vpures=vpure(p,T,**par)
par["vpure"]=vpures
lngi_fun=lambda wi :lngi(T,wi,**par)
dlnai_dlnwi_fun=lambda wi: dlnai_dlnxi(T,wi,**par)

# %% [markdown]
# For a the diffusion of three components, three binary diffusion coefficients need to be defined

# %% [markdown]
# 
# $\hat{Ð} =$ 
# $\left[\begin{array}{rrr} 
# 0 & Ð_{12} & Ð_{13} \\ 
# 0 & 0 & Ð_{23} \\ 
# 0 & 0 & 0 \\ 
# \end{array}\right]$
# 
# $Ð_{vec} = \left[\begin{array}{rrr} Ð_{12} & Ð_{13} & Ð_{23} \end{array}\right]$

# %%
Dvec=np.asarray([1E-13,1E-13,2E-13])
# Dvec=np.asarray([1E-13,1E-13,1E-13])
# Dvec=np.asarray([1E-7,6E-8,1E-14])
# Dvec=np.asarray([8E-9,2E-7,3.5E-10]) #dl 03
# Dvec=np.asarray([8E-10,1E-11,8E-11]) #dl 03
Dvec=np.asarray([8E-10,1E-8,8E-8]) #dl 03

# Dvec=np.asarray([8E-8,2E-7,3.5E-10]) #dl 03 Problem
# Dvec=np.asarray([12E-8,9E-10,11E-12]) #dl 05
# Dvec=np.asarray([8E-5,2E-,3.5E-5]) #dl 05
# Dvec=np.asarray([20E-9,2E-7,3.5E-14])
# Dvec=np.asarray([1E-7,6E-8,1E-14])

# %% [markdown]
# Next we define the time array and which component is mobile

# %%
nt=300
# t=np.linspace(0,16.8,nt)*60
t=np.linspace(0,500,nt)*60
mobile=np.asarray([False,True,True])
# mobile=np.asarray([False,False,True])

# %%
EJ=np.asarray([2E9])
etaJ=np.asarray([1E13])
exponent=np.asarray([0.,2.])
wtid=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,nz=20)
# wt,wtz,zvec,Lt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,full_output=True,nz=20)
wt,wtz,zvec,Lt=Diffusion_MS(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,full_output=True,nz=20,EJ=EJ,etaJ=etaJ,exponent=exponent)
# wt=Diffusion_MS_iter(t,L,Dvec,wi0,wi8,Mi,mobile,swelling=True,dlnai_dlnwi_fun=dlnai_dlnwi_fun)
# notreleased=wt/wi0*Lt[:,None]/L
# release=1-notreleased
# plt.plot((wt*Lt[:,None])[:,1])

wrel=-mges0/msol*(wt*Lt[:,None]/L-wi0)
# wt=m0/msol*release
# wt[:,2]=1-wt[:,0]-wt[:,1]
# wtid=wt


# %% [markdown]
# We can determine the mass dissolved in the dissolution medium by quantifying the mass that leaves the ASD. The initial mass of the ASD and the mass of the dissolution medium must be known

# %% [markdown]
# We plot the results

# %%

font = {'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)

color1 = "#99CC00"
color2 = "#F08208"
color3 = "#99CDE9"

fig, ax = plt.subplots(figsize=(5, 4), dpi = 200)
fig.subplots_adjust(hspace=0.5, wspace=0.3)

# ax.plot(t/60,wtid[:,0], "--",color = color1 , 
#         linewidth = 2.0, label = "wPol")
# ax.plot(t/60,wtid[:,1], "--",color = color2 , 
#         linewidth = 2.0, label = "wAPI")

ax.plot(t/60,wrel[:,0], "-",color = color1 , 
        linewidth = 2.0, label = "wPol")
ax.plot(t/60,wrel[:,1], "-",color = color2 , 
        linewidth = 2.0, label = "wAPI")
ax.plot(t/60,wrel[:,2], "-",color = color3 , 
        linewidth = 2.0, label = "ww")


fig0, ax0 = plt.subplots(figsize=(5, 4), dpi = 200)
ax0.plot(t/60,(Lt-L)*1000, "g-",color = color2 , 
        linewidth = 2.0, label = "wAPI")

ax0.plot(texp,thickness, "gx",color = color2 , 
        linewidth = 2.0, label = "wAPI")

ax.legend(fontsize="small")
ax.set_xlabel('$t$ / s')
ax.set_ylabel('$wi$ / -')
ax.axis([0, 300., 0., 1E-4])
start, end = ax.get_xlim()
# ax.xaxis.set_ticks(np.linspace(start, end, 5))
start, end = ax.get_ylim()
# ax.yaxis.set_ticks(np.linspace(start, end, 5))
plt.show()



# ax2.plot(zexp2*Lt[0],tolu, "o")


# %% [markdown]
# 

# %%
color1 = "#99CC00"
color2 = "#F08208"
color3 = "#99CDE9"
color4 = "#000000"



fig1, ax1 = plt.subplots(figsize=(5, 4), dpi = 200)
fig2, ax2 = plt.subplots(figsize=(5, 4), dpi = 200)
fig3, ax3 = plt.subplots(figsize=(5, 4), dpi = 200)
fig4, ax4 = plt.subplots(figsize=(5, 4), dpi = 200)
fig1.subplots_adjust(hspace=0.5, wspace=0.3)
fig2.subplots_adjust(hspace=0.5, wspace=0.3)
fig3.subplots_adjust(hspace=0.5, wspace=0.3)
fig4.subplots_adjust(hspace=0.5, wspace=0.3)
ax1.set_ylabel('$z$ / µm')
ax1.set_xlabel('$t$ / min')
ax2.set_ylabel('$z$ / µm')
ax2.set_xlabel('$t$ / min')
ax3.set_ylabel('$z$ / µm')
ax3.set_xlabel('$t$ / min')
ax4.set_ylabel('$z$ / µm')
ax4.set_xlabel('$t$ / min')
Micmet,Minutes=np.meshgrid(zvec*1E6,t/60)
expansion=Lt[:,None]/L
# X,Y=Micmet*np.cos(Phi)*expansion,Micmet*np.sin(Phi)*expansion
pl1=ax1.contourf(Minutes,Micmet*expansion,wtz[:,0,:],cmap="Greens")
cbar1 = fig.colorbar(pl1)
pl2=ax2.contourf(Minutes,Micmet*expansion,wtz[:,1,:],cmap="Oranges")
cbar2 = fig.colorbar(pl2)
pl3=ax3.contourf(Minutes,Micmet*expansion,wtz[:,2,:],cmap="Blues")
cbar3 = fig.colorbar(pl3)
pl4=ax4.contourf(Minutes,Micmet*expansion,wtz[:,1,:]/(wtz[:,1,:]+wtz[:,0,:]),cmap="Reds")
cbar4 = fig.colorbar(pl4)

fig12,ax12=plt.subplots()
[ax12.plot(zvec*1E6*Lt[i]/L,wtz[i,0,:], "-",color = color1 , linewidth = 2.0) for i,val in enumerate(wtz[:,0,0])]

fig10 = plt.figure()
ax10 = fig10.add_subplot(141, polar=True)
ax11 = fig10.add_subplot(142, polar=True)
ax12 = fig10.add_subplot(143, polar=True)

phi=np.linspace(0,2*np.pi,41)
Rad,Phi=np.meshgrid(zvec*1E6,phi)
pl10=ax10.contourf(Phi,Rad*expansion[0],np.meshgrid(wtz[0,2,:],phi)[0],cmap="Blues",vmin=0, vmax=1)
# cbar10 = fig.colorbar(pl10)
pl11=ax11.contourf(Phi,Rad*expansion[100],np.meshgrid(wtz[100,2,:],phi)[0],cmap="Blues",vmin=0, vmax=1)
# cbar11 = fig.colorbar(pl11)
pl12=ax12.contourf(Phi,Rad*expansion[-1],np.meshgrid(wtz[-1,2,:],phi)[0],cmap="Blues",vmin=0, vmax=1)
# cbar12 = fig.colorbar(pl12)

ax10.grid(False)
ax11.grid(False)
ax12.grid(False)

ax10.set_xticklabels([])
ax10.set_yticklabels([])
ax11.set_xticklabels([])
ax11.set_yticklabels([])
ax12.set_xticklabels([])
ax12.set_yticklabels([])

ax10.set_ylim(0, np.max(zvec*1E6*expansion))
ax11.set_ylim(0, np.max(zvec*1E6*expansion))
ax12.set_ylim(0, np.max(zvec*1E6*expansion))

ax10.spines['polar'].set_visible(False)
ax11.spines['polar'].set_visible(False)
ax12.spines['polar'].set_visible(False)

ax30 = fig10.add_subplot(1,4,4)
# ax30.axis('off')
sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin=0, vmax=1))
plt.colorbar(sm,ax30)


circular(t,zvec,wtz,Lt=Lt,instances=6,comp=2,cmap="Blues")
# cbar = plt.colorbar(pl12,ax30)

# ax1.plot(zexp1*Lt[0],meth, "o")

# [ax2.plot(zvec*1E6*Lt[i]/L,wtz[i,1,:], "-",color = color2 , linewidth = 2.0) for i,val in enumerate(wtz[:,0,0])]




# [ax3.plot(zvec*1E6*Lt[i]/L,wtz[i,2,:], "-",color = color3 , linewidth = 2.0) for i,val in enumerate(wtz[:,0,0])]

# [ax4.plot(zvec*1E6*Lt[i]/L,wtz[i,1,:]/(1-wtz[i,2,:]), "-",color = color3 , linewidth = 2.0) for i,val in enumerate(wtz[:,0,0])]


# [ax3.plot(zvec*1E6,wtz[i,1,:]/(wtz[i,1,:]+wtz[i,0,:]), "-",color = color3 , linewidth = 2.0) for i,val in enumerate(wtz[:,0,0])]
plt.show()