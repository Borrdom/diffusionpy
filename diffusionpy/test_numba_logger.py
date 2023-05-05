from numba import njit
import numpy as np

def logging_jit(func):
    def inner(*args, **kwargs):
        origsigs = set(func.signatures)
        result = func(*args, **kwargs)
        newsigs = set(func.signatures)
        if newsigs != origsigs:
            new = (newsigs ^ origsigs).pop()
             # PUT YOUR LOGGER HERE!
            print("Numba compiled for signature: {}".format(new))
        return result
    return inner

@logging_jit
@njit
def foo(a):
    return a + 1

print(foo(4)) # first compile and run for int
print(foo(5)) # second is just a run, int sig is cached
print(foo(6.7)) # third is a compile and run for float