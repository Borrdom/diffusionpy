import numpy as np
x=np.linspace(0,10,100)
y=x**2
dy1=np.gradient(y,edge_order=2)
dy2=np.convolve(y,[1,-1])

