# diffusionpy

![alt text](https://github.com/Borrdom/diffusionpy/blob/main/PyFusion.png?raw=true)

Provides a one-dimensional Stefan-Maxwell diffusion model with a PC-SAFT implementation

# Installation
```console
git clone https://github.com/Borrdom/diffusionpy
cd diffusionpy
pip install .
```
# Minimal working example

If no errors occoured during installation, the following lines will model the diffusion kinetics of three substances from a starting mass fraction of wi0 to the mass fraction wi8
```python
import numpy as np
from diffusionpy import Diffusion_MS
import matplotlib.pyplot as plt
t=np.linspace(0,300,30) # time points in seconds
D=np.asarray([1E-12,1E-12,5E-14]) # diffusion coefficients in meters per second
wi0=np.asarray([0.33,0.33,0.33]) # mass fractions at t=0
wi8=np.asarray([0.01,0.01,0.98]) # mass fractions in equilibrium
L=5E-6 # diffusion path in meters
mobile=np.asarray([True,True,False]) # specify mobile components
wt=Diffusion_MS(t,L,D,wi0,wi8,mobile)[0]
fig,ax=plt.subplots()
ax.plot(t,wt)
ax.set_xlabel("time / s") 
ax.set_ylabel("mass fractions / -") 
ax.set_xticks(np.linspace(0,300,5))
ax.set_yticks(np.linspace(0,1,5))
plt.show()
```

# Getting started
check out the jupyter notebooks in the [examples](https://github.com/Borrdom/diffusionpy/tree/main/examples_notebooks) folder.


# Documentation
can be found [here](https://github.com/Borrdom/diffusionpy/tree/main/docs/html/index.html).

# License information

BSD-3

# Originlab-like plot format 

Just copy the https://github.com/Borrdom/diffusionpy/tree/main/.matplotlib folder into your %USERPROFILE%. This will alter the matplotlib rc permantly so that the format below will be your new default. 

```python
import matplotlib.pyplot as plt
import numpy as np
x=np.linspace(1,10,11)
y1=x+1
y2=x+2
y3=x+3
fig,ax=plt.subplots()
ax.plot(x,x,'ko-')
ax.plot(x,y1,'C6s-')
ax.plot(x,y2,'C3^-')
ax.plot(x,y3,'C0*-')
ax.set_xlabel(r'$m_{Aminoacid}/(gmol^{-1})$')
ax.set_ylabel(r'$\alpha/(Jmol^{-1})$')
ax.set_xticks(np.linspace(0,12,7))
ax.set_yticks(np.linspace(0,14,8))
plt.show()
```
![alt text](https://github.com/Borrdom/diffusionpy/blob/main/originlike.png?raw=true)

Moreover, Matplotlib will recognize the following formatstrings as these colors.


- C0 : ![#ff1778](https://via.placeholder.com/15/99CC00/000000?text=+) `#99CC00`
- C1 : ![#00dbb2](https://via.placeholder.com/15/99CDE9/000000?text=+) `#99CDE9`
- C2 : ![#69adff](https://via.placeholder.com/15/246FE2/000000?text=+) `#246FE2`
- C3 : ![#9ece6a](https://via.placeholder.com/15/FF8500/000000?text=+) `#FF8500`
- C4 : ![#c7ccd4](https://via.placeholder.com/15/FFCCCC/000000?text=+) `#FFCCCC`
- C5 : ![#666b88](https://via.placeholder.com/15/FFD67E/000000?text=+) `#FFD67E`
- C6 : ![#0c0d12](https://via.placeholder.com/15/666666/000000?text=+) `#666666`

