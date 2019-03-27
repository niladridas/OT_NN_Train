#!/usr/bin/python
from dualDyn import param_dyn, meas, state_dyn
import numpy as np


w_curr = np.array([10,8/3,28])
v_curr = 0
mod_state_dyn = lambda x_prev: state_dyn(x_prev,w_curr,v_curr) # Using lambda function
x_prev = np.array([1,1,1])
a = mod_state_dyn(x_prev)
print a
