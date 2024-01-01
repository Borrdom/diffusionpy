import autograd.numpy as np  # Thinly-wrapped numpy
from autograd import elementwise_grad    # The only autograd function you may ever need

def tanh(x):                 # Define a function
    y = np.exp(-2.0 * x)
    return (1.0 - y) / (1.0 + y)

grad_tanh = elementwise_grad(tanh)       # Obtain its gradient function
grad_tanh(1.0)               # Evaluate the gradient at x = 1.0

