import numpy as np
days=22
ran=np.random.rand(days)
ub=9
lb=7
delta=ub-lb
rand=ran*delta+lb
hou=np.arange(lb,ub+0.00001,0.5)
final=np.zeros(days)
for i in range(days):
    delt=abs(rand[i]-hou)
    ind=np.argmin(delt)
    final[i]=hou[ind]
print(np.mean(final))    