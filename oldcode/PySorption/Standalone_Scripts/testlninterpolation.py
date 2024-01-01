import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
x=np.linspace(0,1,5)
x2=np.linspace(0,1,50)

y=8**x
y2=8**x2

y_fun=InterpolatedUnivariateSpline(x,y,k=1)
ysim=y_fun(x2)
ylog_fun=InterpolatedUnivariateSpline(x,np.log(y),k=1)
ysimlog=np.exp(ylog_fun(x2))
#plt.plot(x,y,'kx')
plt.plot(x2,ysimlog,'r--')
plt.plot(x2,ysim,'k-.')

plt.plot(x2,y2,'kx')
plt.show()
