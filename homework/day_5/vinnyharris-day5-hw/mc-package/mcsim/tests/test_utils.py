"""
Test for mcsim package - utils module.
"""
import math
import mcsim.utils as utils
import random

def test_accept_or_reject_deltaneg():
    """
    Test calculate distance function
    """
    delta = -1
    beta = 1

    observed = utils.accept_or_reject(delta, beta)
    
    assert observed == True

def test_accept_or_reject_random1():
    """
    Test calculate distance function
    """
    random.seed(0)
    delta = 1
    beta = 1

    observed = utils.accept_or_reject(delta, beta)
    
    assert observed == True

def test_accept_or_reject_random1():
    """
    Test calculate distance function
    """
    random.seed(1)
    delta = 1
    beta = 1

    observed = utils.accept_or_reject(delta, beta)
    
    assert observed == True

def test_tail_correction():
    """
    Test calculate tail correction
    """
    expected = -198.4888837441566
    observed = utils.calculate_tail_correction(10.0, 800, 3.0)

    assert math.isclose(expected, observed)
