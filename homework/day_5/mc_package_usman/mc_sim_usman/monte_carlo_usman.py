"""
Functions for running a Motne Carlo Simulation
"""

import math
import random
import mc_sim_usman.utilities_usman as ut

def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Calculate the total energy of a set of particles using the Lennard Jones potential.
    
    Parameters
    ----------
    coordinates : list
        A nested list containing the x, y,z coordinate for each particle
    box_length : float
        The length of the box. Assumes cubic box.
    cutoff : float
        The cutoff length
    
    Returns
    -------
    total_energy : float
        The total energy of the set of coordinates.
    """
    
    total_energy = 0
    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            # Calculate the distance between the particles - exercise.
            dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length)

            if dist_ij < cutoff:
                # Calculate the pairwise LJ energy
                LJ_ij = calculate_LJ(dist_ij)

                # Add to total energy.
                total_energy += LJ_ij
    return total_energy


def calculate_LJ(r_ij):
    """
    The LJ interaction energy between two particles.

    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.
    
    Returns
    -------
    pairwise_energy : float
        The pairwise Lennard Jones interaction energy in reduced units.

    Examples
    --------
    >>> calculate_LJ(1)
    0

    """
    
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    
    return pairwise_energy


def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two points. When box_length is set, the minimum image convention is used to calculate the distance between the points.

    Parameters
    ----------
    coord1, coord2 : list
        The coordinates of the points, [x, y, z]
    
    box_length : float, optional
        The box length

    Returns
    -------
    distance : float
        The distance between the two points accounting for periodic boundaries
    """
    distance = 0
        
    for i in range(3):
        hold_dist = abs(coord2[i] - coord1[i])
    
        if (box_length):    
            if hold_dist > box_length/2:
                hold_dist = hold_dist - (box_length * round(hold_dist/box_length))
        distance += math.pow(hold_dist, 2)

    return math.sqrt(distance)


def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system)
    
    Parameters
    ----------
    coordinates : list
        The coordinates for all the particles in the system.
        
    i_particle : int
        The particle number for which to calculate the energy.
        
    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated.
    
    box_length : float
        The length of the box for periodic bounds
        
    Returns
    -------
    e_total : float
        The pairwise interaction energy of the ith particles with all other particles in the system
    """
    
    e_total = 0.0
    #creates a list of the coordinates for the i_particle
    i_position = coordinates[i_particle]
    
    num_atoms = len(coordinates)
    
    for j_particle in range(num_atoms):
        
        if i_particle != j_particle:
            #creates a list of coordinates for the j_particle
            j_position = coordinates[j_particle]
            rij = calculate_distance(i_position, j_position, box_length)
            
            if rij < cutoff:
                e_pair = calculate_LJ(rij)
                e_total += e_pair
    
    return e_total
