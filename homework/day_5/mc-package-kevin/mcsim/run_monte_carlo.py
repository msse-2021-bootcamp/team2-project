"""
Runs the Monte Carlo Simulation using the directed implementation
"""
import math
import random
import monte_carlo as mc
import monte_carlo_numpy as numpy
def run_simulation(run_simulation_np(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000, implement=None):
    if implement == None:
        mc.run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000)
    elif implement == numpy:
        implement.numpy.run_simulation_np(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000)
    else:
        print("Implementation not recognized in this version.")
