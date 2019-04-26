# curr: time step k
# prev: time step k-1
# next: time step k+1
# Lorenz System w=[sigma; rho; beta] x=[x;y;z]
import numpy as np

def param_dyn(w_prev, u_curr): # Eqn.21
    w_curr = w_prev + u_curr
    return w_curr;

def state_dyn(x_prev,w_curr,v_curr): # Eqn.19
    dt = 1 # Time discretization interval
    xk = x_prev[0] + dt*w_curr[0]*(x_prev[1]-x_prev[0])
    yk = x_prev[1] + dt*(x_prev[0]*(w_curr[1] - x_prev[2]) - x_prev[1])
    zk = x_prev[2] + dt*(x_prev[0]*x_prev[1] - w_curr[2]*x_prev[2])
    x_curr = np.array([xk, yk, zk])
    return x_curr;

def meas(x_curr,n_curr): # Eqn.22
    y_curr = x_curr + n_curr
    return y_curr;
