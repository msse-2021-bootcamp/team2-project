"""
Test for mcsim package - monte carlo module.
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


# write a test for the calculate distance function which tests for the periodic boundary conditions
def test_calculate_distance_2():
    """
    Test periodic boundary condition for calculate distance function 
    """
    point1=[0,0,0]
    point2=[0,0,8]
    box_length = 10 
    expected = 2
    observed = mc.calculate_distance(point1,point2,box_length)
    assert expected == observed


def test_calculate_LJ_1():
    """
    Test the Lennard Jones pair energy 
    """
    assert mc.calculate_LJ(1) == 0


def test_calculate_LJ_2():
    """
    Test the Lennard Jones pair energy 
    """
    assert mc.calculate_LJ(math.pow(2, (1/6))) == -1.0


def test_calculate_total_energy():
    """
    Test calculating the total energy 
    """
    coordinates = [[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]]

    assert mc.calculate_total_energy(coordinates, 10, 3.0) == -2.031005859375


def test_calculate_total_energy_1():
    """
    Test calculating the total energy 
    """
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]

    assert mc.calculate_pair_energy(coordinates, 1, 10, 3) == -2

def test_calculate_total_energy_2():
    """
    Test calculating the total energy 
    """
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]

    assert mc.calculate_pair_energy(coordinates, 0, 10, 3) == mc.calculate_pair_energy(coordinates, 2, 10, 3)

def test_calculate_total_energy_3():
    """
    Test calculating the total energy 
    """
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]
    assert mc.calculate_pair_energy(coordinates, 0, 10, 2) == -1
