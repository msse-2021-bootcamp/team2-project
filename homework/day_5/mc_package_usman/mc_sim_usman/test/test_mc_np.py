"""
Test for mc_sim_usman package.
"""
import math
import random
import numpy as np
import mc_sim_usman as mcn


def test_calculate_distance_1():
    """
    Test calculate distance function
    """
    point1=[0,0,0]
    point2=[0,1,0]

    expected = 1
    observed = mcn.calculate_distance(point1,point2)
    assert math.isclose(expected, observed)

# write a test for the calculate distnace funtion which tests for the periodic boundary conditions

def test_calculate_distance_2():
    """
    Test periodic boundary condition for calculate distance function 
    """
    point1=[0,0,0]
    point2=[0,0,8]
    box_length = 10 
    expected =2
    observed = mcn.calculate_distance(point1,point2,box_length)
    assert expected == observed

def test_calculate_distance_np_1():
    """
    Test calculate distance function
    """
    point1=np.array([0,0,0])
    point2=np.array([0,1,0])

    expected = 1
    observed = mcn.calculate_distance_np(point1,point2)
    assert math.isclose(expected, observed)

# write a test for the calculate distnace funtion which tests for the periodic boundary conditions

def test_calculate_distance_np_2():
    """
    Test periodic boundary condition for calculate distance function 
    """
    point1=np.array([0,0,0])
    point2=np.array([0,0,8])
    box_length = 10 
    expected =2
    observed = mcn.calculate_distance_np(point1,point2,box_length)
    assert expected == observed

def test_calculate_LJ():
    """
    Test calculate LJ function
    """
    r_ij = 1
    expected = 0
    observed = mcn.calculate_LJ(r_ij)
    assert math.isclose(expected, observed)
    r_ij = math.pow(2,(1/6))
    expected = -1
    observed = mcn.calculate_LJ(r_ij)
    assert math.isclose(expected, observed)

# write a test for the calculate distnace funtion which tests for the periodic boundary conditions

def test_calculate_LJ_np():
    """
    Test calculate LJ function
    """
    r_ij = np.array([1])
    expected = 0
    observed = mcn.calculate_LJ_np(r_ij)
    assert math.isclose(expected, observed)

    r_ij = np.array([math.pow(2,(1/6))])
    expected = -1
    observed = mcn.calculate_LJ_np(r_ij)
    assert math.isclose(expected, observed)

    r_ij = np.array([1, 2])
    expected = -0.0615234375
    observed = mcn.calculate_LJ_np(r_ij)
    assert math.isclose(expected, observed)

    r_ij = np.array([1, 4, 3])
    expected = -0.006455765825659674
    observed = mcn.calculate_LJ_np(r_ij)
    assert math.isclose(expected, observed)

def test_calculate_total_energy():
    """
    Test calculate total energy function
    """
    coordinates = [[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]]

    expected = -2.031005859375

    observed = mcn.calculate_total_energy(coordinates, 10, 3.0) 
    
    assert math.isclose(expected, observed)

# write a test for the calculate distnace funtion which tests for the periodic boundary conditions

def test_calculate_total_energy_np():
    """
    Test calculate total energy np function
    """
    coordinates = [[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]]

    expected = -2.031005859375

    observed = mcn.calculate_total_energy_np(coordinates, 10, 3.0) 
    
    assert math.isclose(expected, observed)

def test_calculate_pair_energy():
    """
    Test calculate total energy function
    """
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]]

    expected = -2
    observed = mcn.calculate_pair_energy(coordinates, 1, 10, 3) 
    assert math.isclose(expected, observed)

    expected = -1
    observed = mcn.calculate_pair_energy(coordinates, 0, 10, 2) 
    assert math.isclose(expected, observed)


# write a test for the calculate distnace funtion which tests for the periodic boundary conditions

def test_calculate_pair_energy_np():
    """
    Test calculate pair energy np function
    """
    coordinates = [[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]]
    
    expected = -2
    observed = mcn.calculate_pair_energy_np(coordinates, 1, 10, 3.0) 
    assert math.isclose(expected, observed)

    expected = -1
    observed = mcn.calculate_pair_energy_np(coordinates, 0, 10, 2) 
    assert math.isclose(expected, observed)

def test_calculate_test_correction():
    """
    Test calculate test correction function
    """
    num_of_particles = 800
    box_length = 10
    cutoff = 3

    expected = -198.4888837441566
    observed = mcn.calculate_tail_correction(num_of_particles, box_length, cutoff)
    assert math.isclose(expected, observed)

def test_accept_or_reject():
    """
    Test to calculate the accept or reject function
    """
    beta = 1
    delta_U = -5
    expected = True
    observed = mcn.accept_or_reject(delta_U, beta)
    assert observed == expected
    
    beta = 1
    delta_U = 0
    expected = True
    observed = mcn.accept_or_reject(delta_U, beta)
    assert observed == expected




