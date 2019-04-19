# Unscented Kalman Filter
UKF in Python for **State, Parameter, and Dual Estimation**
****
This work is based on:
**[The Unscented Kalman Filter for Nonlinear Estimation]( https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf)**
****
Contributors:
1. Vaishnav Tadiparthi
2. Niladri Das
3. Vedang Deshpande
****
**Dual Estimation:**<br />
Discrete-time Nonlinear Dynamics for dual estimation problem
$$\mathbf{x}_{k+1} = \mathbf{F}(\mathbf{x}_k,\mathbf{v}_k,\mathbf{w})\\
\mathbf{y}_k = \mathbf{H}(\mathbf{x}_k,\mathbf{n}_k,\mathbf{w})$$
where the system states are denoted by $\mathbf{x}_k$.
