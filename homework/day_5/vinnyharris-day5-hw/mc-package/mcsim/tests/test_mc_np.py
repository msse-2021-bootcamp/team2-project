"""
Test for mcsim package.
"""
import math
import mcsim.monte_carlo_np as mc
import numpy as np


def test_calculate_distance_np_1():
    """
    Test calculate distance function
    """
    point1=np.array([0,0,0])
    point2=np.array([0,1,0])

    expected = np.array([1.0])
    observed = mc.calculate_distance_np(point1,point2)
    assert np.array_equal(expected, observed)


# write a test for the calculate distance function which tests for the periodic boundary conditions
def test_calculate_distance_np2():
    """
    Test periodic boundary condition for calculate distance function 
    """
    point1= np.array([[0, 0, 0],[0, 1, 0]])
    point2= np.array([[0, 8, 0],[0, 1.5, 0]])
   
    box_length = 3.3
    expected = np.array([1.4, 0.5])
    observed = mc.calculate_distance_np(point1,point2,box_length)
    assert np.allclose(expected, observed)


def test_calculate_lj_np_1():
    """
    Test the Lennard Jones pair energy 
    """
    assert mc.calculate_lj_np(1) == 0


def test_calculate_lj_np_2():
    """
    Test the Lennard Jones pair energy 
    """
    assert mc.calculate_lj_np(math.pow(2, (1/6))) == -1.0


def test_calculate_total_energy_np():
    """
    Test calculating the total energy 
    """
    coordinates = np.array([[0, 0, 0], [0, math.pow(2, 1/6), 0], [0, 2*math.pow(2, 1/6), 0]])

    assert np.isclose(mc.calculate_total_energy_np(coordinates, 10, 3.0), -2.031005859375)


def test_calculate_pair_energy_np_1():
    """
    Test calculating the total energy 
    """
    coordinates = np.array([[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]])

    assert mc.calculate_pair_energy_np(coordinates, 1, 10, 3) == -2

def test_calculate_pair_energy_np_2():
    """
    Test calculating the total energy 
    """
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]
    assert mc.calculate_pair_energy_np(coordinates, 0, 10, 3) == mc.calculate_pair_energy_np(coordinates, 2, 10, 3)

def test_calculate_pair_energy_3():
    """
    Test calculating the total energy 
    """
    coordinates = [[0, 0, 0], [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]
    assert mc.calculate_pair_energy_np(coordinates, 0, 10, 2) == -1
