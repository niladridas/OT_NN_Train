#!/usr/bin/python


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
def weightpropdualukf(postsamples,weights):
  """ Weight update function for Dual UKF

  :return: propagated prior samples, expected observations
  """
  

def weightupdatedualukf():

# TO-DO: Signal filter

