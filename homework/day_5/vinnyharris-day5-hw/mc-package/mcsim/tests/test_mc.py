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

def test_calculate_distnace_2():
    """
    Test periodic boundary condition for calculate distance function 
    """
    point1=[0,0,0]
    point2=[0,0,8]
    box_length = 10 
    expected =2
    observed = mc.calculate_distance(point1,point2,box_length)
    assert expected == observed