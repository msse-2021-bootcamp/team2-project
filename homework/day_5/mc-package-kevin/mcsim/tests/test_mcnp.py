"""
Test for mcsim package, monte_carlo_numpy module.
"""
import math
import numpy as np
import mcsim.monte_carlo_numpy as mcnp

def test_calculate_distance_1():
    """
    Test calculate distance function
    """
    test_coords = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0,3]])
    expected = np.array([0,1,2,3])
    observed = mcnp.calculate_distance_np(test_coords[0],test_coords)
    assert np.array_equal(expected,observed)

# write a test for the calculate distnace funtion which tests for the periodic boundary conditions

def test_calculate_distance_2():
    """
    Test periodic boundary condition for calculate distance function 
    """
    test_coords = np.array([[0, 0, 0], [3, 0, 0], [6, 0, 0], [0, 0, 9]])
    box_length = 10 
    expected = np.array([0,3,4,1])
    observed = mcnp.calculate_distance_np(test_coords[0],test_coords,box_length)
    assert np.array_equal(expected,observed)

def test_calculate_LJ():
    """ 
    Test calculate LJ function:
    """
    test_coords = np.array([1,2**(1/6)])
    expected = np.array([0,-1])
    observed = mcnp.calculate_LJ_np(test_coords)
    assert np.array_equal(expected,observed)

def test_calculate_pair_energy_1():
    """
    Test calculate pair energy function.
    """
    test_coords = np.array([[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]])
    box_length, cutoff=10,3
    expected= -2
    observed= mcnp.calculate_pair_energy_np(test_coords, 1,  box_length, cutoff)
    assert expected == observed

def test_calculate_pair_energy_2():
    """
    Test cutoff condition for calculate pair energy function.
    """
    test_coords = np.array([[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]])
    box_length, cutoff=10,2
    expected= -1
    observed= mcnp.calculate_pair_energy_np(test_coords, 0, box_length, cutoff)
    assert expected == observed

def test_calculate_total_energy():
    """
    Test the calculate total energy function
    """
    test_coords = np.array([[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]])
    box_length, cutoff=10,2
    expected = -2
    observed= mcnp.calculate_total_energy_np(test_coords, box_length, cutoff)
    assert expected == observed

