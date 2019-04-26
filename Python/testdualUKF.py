#!/usr/bin/python
import numpy as np


def main():
    """ Dual UKF for Discrete time systems

    This function test the dual UKF algorithm.
    Two UKF algorithms are run simultaneously.
    The first algorithm is for updating the weights.
    These updates weights are used in the next UKF algorithm to update the states.
    The updated states are then fed back to the first UKF algorithm.
    """
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
  main()


#TO-DO: Observation sequence
# yobs = # An array with dimension no_of_y_states X no_of_observations

# TO-DO: Weight filter
def weightpropdualukf(wsigmapts, sigw, mu, xest, xprocesscov, ymeascov):
  """ Weight update function for Dual UKF

  : params: mu : puu = mu*pww (from paper), for data length of 1000, mu = 1e-4
  :return: propagated prior samples, expected observations
  """
  # The weight sigma points are propagate using the equation
  # W_next = w_previous + u
  # u: zero mean noise with cov Pu = mu*Pw
  # TO-DO: Pw: calculated from the wsigmapts (not clear from the paper)
  l1 = np.shape(wsigmapts)[0]
  wm = np.zeros((l1, 1))
  for i in range(0, 2 * l1):
    wm = wm + sigw.wm[i]*wsigmapts[:, i]
  for i in range(0, 2 * l1):
    pww = pww + sigw.wc[i] * np.matmul(wsigmapts[:, i] - wm, (wsigmapts[:, i] - wm).transpose())
  puu = mu*pww
  meanu = np.zeros((1, l1))
  covu = puu
  for i in range(0, 2 * l1):
    u = np.random.multivariate_normal(meanu, covu, (3, 3))
    wsigmapts[:, i] = wsigmapts[:, i] + u
  for i in range(0, 2 * l1):


def weightupdatedualukf():

# TO-DO: Signal filter

