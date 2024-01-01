import matplotlib.pyplot as plt
import numpy as np

T=298.15
Tg=np.linspace(288.15,348.15)

eta0=0.36*5E9
logeta1=np.log10(8E12/eta0)
logeta2=np.log10(3E11/eta0)
Tg1=48-298.15
Tg2=-10-298.15 #27

C1=14
C2=105
C3=(logeta2-logeta1)/(Tg2/T-Tg1/T)
C4=logeta2+C3*Tg2/T
WL=-C1*(T-Tg)/(C2+T-Tg)
AR=C3*Tg/T+C4
fig,ax=plt.subplots()
ax.plot(T-Tg,WL+np.log10(eta0))
ax.plot(T-Tg,AR+np.log10(eta0))
plt.show()
