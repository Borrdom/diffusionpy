import numpy as np

def weno_interpolation(u):
    n = len(u)
    u_interp = np.zeros_like(u)
    
    # WENO parameters
    stencil_size = 3  # Size of the stencil
    smoothness_eps = 1e-6  # Small constant to avoid division by zero
    
    for i in range(n):
        # Select the stencil
        stencil_start = max(i - (stencil_size - 1) // 2, 0)
        stencil_end = min(stencil_start + stencil_size, n)
        stencil = u[stencil_start:stencil_end]
        
        # Calculate smoothness indicators
        smoothness = np.zeros_like(stencil)
        smoothness[1:-1] = (1/6) * (stencil[2:] - 2*stencil[1:-1] + stencil[:-2])**2 + \
                           (2/3) * (stencil[2:] - stencil[:-2])**2 + \
                           smoothness_eps
        
        # Calculate weights
        weights = np.zeros_like(stencil)
        weights[1:-1] = (1/10) / ((smoothness[1:-1] + smoothness_eps)**2)
        weights /= np.sum(weights)
        
        # Interpolate the solution
        u_interp[i] = np.sum(weights * stencil)
    
    return u_interp

# Example usage
u = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
u_interp = weno_interpolation(u)
print(u_interp)