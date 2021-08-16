"""
Test for mcsim package.
"""
import math
import mcsim.monte_carlo as mc

def test_calculate_distance_1():
    """
    Test calculate distance function
    """
    point1=[0,0,0]
    point2=[0,1,0]

    expected = 1
    observed = mc.calculate_distance(point1,point2)
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
    observed = mc.calculate_distance(point1,point2,box_length)
    assert expected == observed

def test_calculate_LJ1():
    """ 
    Test calculate LJ function:
    """
    distance=1
    expected = 0
    observed = mc.calculate_LJ(distance)
    assert expected == observed

def test_calculate_LJ2():
    """ 
    Test calculate LJ function:
    """
    distance = 2**(1/6)
    expected = -1
    observed = mc.calculate_LJ(distance)
    assert expected == observed

def test_calculate_pair_energy_1():
    """
    Test calculate pair energy function.
    """
    test_coords = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]]
    box_length, cutoff=10,3
    expected= -2
    observed= mc.calculate_pair_energy(test_coords, 1,  box_length, cutoff)
    assert expected == observed

def test_calculate_pair_energy_2():
    """
    Test cutoff condition for calculate pair energy function.
    """
    test_coords = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]]
    box_length, cutoff=10,2
    expected= -1
    observed= mc.calculate_pair_energy(test_coords, 0, box_length, cutoff)
    assert expected == observed

def test_calculate_total_energy():
    """
    Test the calculate total energy function
    """
    test_coords = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*2**(1/6)]]
    box_length, cutoff=10,2
    expected = -2
    observed= mc.calculate_total_energy(test_coords, box_length, cutoff)
    assert expected == observed
