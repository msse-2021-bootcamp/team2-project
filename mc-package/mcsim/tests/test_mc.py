"""
Tests for mcsim package

"""
import math
import mcsim.monte_carlo as mc

def test_calculate_distance():
    point_1 = [0, 0, 0]
    point_2 = [0, 1, 0]

    expected = 1
    observed = mc.calculate_distance(point_1, point_2)
    
    assert math.isclose(expected, observed)


    


# Write a test for the calculate_distance function for periodic boundary conditions
def test_calculate_distance_periodic():
    point_1 = [0, 8, 0]
    point_2 = [0, 0, 0]
    box_length = 10

    expected_distance = 2
    dist1 = mc.calculate_distance(point_1, point_2, box_length)
    assert dist1 == expected_distance