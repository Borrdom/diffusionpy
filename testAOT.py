import numpy as np
from numba.pycc import CC
from numba import njit

cc = CC('my_module')

@njit(['f8[::1](f8[::1],f8[::1])'],cache=True)
def add_JIT(a,b):
    return a+b

@cc.export('add_AOT','f8[::1](f8[::1],f8[::1])')
def add_AOT(a,b):
    return a+b


cc.compile()
import my_module

a=np.asarray([1.,2.])
b=np.asarray([2.,3.])
c1=add_JIT(a,b)
c2=my_module.add_AOT(a,b)


print(c1)
print(c2)