import casadi as cs
import numpy as np
import time
import random
# import matplotlib.pyplot as plt

def Model(input1,input2,input3):
    x=cs.SX.sym("x",2)
    x0=[input1,input2]
    ode=cs.SX.zeros(2)


    ode[1]=x[0]*input3
    ode[0]=-x[1]
    tspan=np.linspace(0,10,100)
    opts={"grid":tspan}
    dae={"x":x,"ode":ode}
    X=cs.integrator("X","cvodes",dae,opts)
    
    if random.random()>0.9:
        time.sleep(2)
    elif random.random()>0.95:
        raise Exception("Some Error")
        
    x_int=cs.horzcat(cs.DM(x0),X(x0=x0)["xf"]).full()

    # fig,ax=plt.subplots()
    # ax.plot(tspan,x_int[0,:])
    # ax.plot(tspan,x_int[1,:])

    return tspan,x_int
