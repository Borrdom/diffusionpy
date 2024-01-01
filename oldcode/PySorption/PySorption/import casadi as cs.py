import casadi as cs
x=cs.SX.sym("x")
y=cs.SX.sym("y")
z=cs.conditional(0,[x],y)
print(z)