# Takes as input 4 functions:
# 1. State dynamics
# 2. State measurement model
# 3. Parameter dynamics
# 4. Parameter measurement model
#
from scipy.integrate import ode


# Ode integration
r = ode(f, jac).set_integrator('zvode', method='bdf', with_jacobian=True)
r.set_initial_value(y0, t0).set_f_params(2.0).set_jac_params(2.0)
t1 = 10
dt = 1
while r.successful() and r.t < t1:
    r.integrate(r.t+dt)
    print("%g %g" % (r.t, r.y))
