import casadi as cs
def matmultadd(a):
    return (a*cs.SX.ones((a.shape[0],a.shape[0]))).T+a
a=cs.SX(3,6)
mi=cs.SX(4,1)
