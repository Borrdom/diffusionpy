from scipy.integrate import solve_ivp
import numpy as np
from matplotlib import pyplot as plt

def bouncing_ball_rhs(t, Y, g=9.81, restitution=0.9):
    """
    Right-hand side of the ODEs for the bouncing ball.
    """
    x, y, vx, vy = Y
    dxdt = vx
    dydt = vy
    dvxdt = 0  # No horizontal acceleration
    dvydt = -g #if y > 0 else -g * restitution  # Change direction and reduce speed on bounce
    return dxdt, dydt, dvxdt, dvydt

def hit_ground(t, Y, g=9.81, restitution=0.9):
    """
    Event function triggered when the ball hits the ground.
    Adjusts the vertical velocity based on the coefficient of restitution.
    """
    _, y, _, vy = Y
    # if y <= 0:
    #     vy *= -restitution
    return y
hit_ground.terminal = True
hit_ground.direction = -1

def main():
    # Initial conditions: [x, y, vx, vy]
    ic = [0, 10, 0, 10]  # Starting at rest
    t_span = (0, 20)  # Simulation time span
    tspan=np.linspace(0,20,50)
    g = 9.81  # Gravity
    restitution = 0.9  # Coefficient of restitution for the bounce

    sol = solve_ivp(bouncing_ball_rhs, t_span, ic,t_eval=tspan, method='LSODA',  events=hit_ground, rtol=1e-6, atol=1e-9,dense_output=True)
    plt.figure(figsize=(8, 6))
    plt.plot(sol.t, sol.y[1])

    t_span = (sol.t[-1],t_span[-1]) 
    tspan=np.linspace(sol.t[-1],t_span[-1],50)
    sol = solve_ivp(bouncing_ball_rhs,t_span, -sol.y[:,-1]*0.8,t_eval=tspan, method='LSODA',  events=hit_ground, rtol=1e-6, atol=1e-9,dense_output=True)
    plt.plot(sol.t, sol.y[1])
    
    t_span = (sol.t[-1],t_span[-1]) 
    tspan=np.linspace(sol.t[-1],t_span[-1],50)
    sol = solve_ivp(bouncing_ball_rhs,t_span, -sol.y[:,-1]*0.8, t_eval=tspan,method='LSODA',  events=hit_ground, rtol=1e-6, atol=1e-9,dense_output=True)
    t_span = (sol.t[-1],t_span[-1]) 
    tspan=np.linspace(sol.t[-1],t_span[-1],50)
    plt.plot(sol.t, sol.y[1])
    sol = solve_ivp(bouncing_ball_rhs,t_span, -sol.y[:,-1]*0.8, t_eval=tspan,method='LSODA',  events=hit_ground, rtol=1e-6, atol=1e-9,dense_output=True)
 
    plt.plot(sol.t, sol.y[1])
    # Plotting

    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.title('Bouncing Ball Simulation')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()