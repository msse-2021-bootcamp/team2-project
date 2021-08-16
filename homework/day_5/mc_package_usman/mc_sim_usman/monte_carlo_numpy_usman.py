"""
Functions for running a Motne Carlo Simulation using NumPy
"""

import numpy as np
import math
import random
import mc_sim_usman.utilities_usman as ut
import time

def calculate_distance_np(coord1, coord2, box_length=None):
    """
    Calculate the distance between two points. When box_length is set, the minimum image convention is used to calculate the distance between the points.

    Parameters
    ----------
    coord1, coord2 : np.array
        The coordinates of the points, [x, y, z]
    
    box_length : float, optional
        The box length

    Returns
    -------
    distance : np.array
        The array of distances between 1 point and the rest accounting for periodic boundaries
    """
    coord_dist=coord1-coord2
    if box_length:    
            coord_dist = coord_dist - (box_length * np.round(coord_dist/box_length))
    if coord_dist.ndim <2:
        #reshaping the array, make it have 2 dimensions
        coord_dist = coord_dist.reshape (1,-1)
    coord_dist=coord_dist**2
    dist_sum= coord_dist.sum(axis=1)
    distance = np.sqrt(dist_sum)
    return distance

def calculate_LJ_np(r_ij):
    """
    The LJ interaction energy between two particles.

    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : np.array()
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
    #we get a list of r_ij from the distance function in the form of an np.array()
    
    r6_term = (1./r_ij) ** 6
    r12_term = (1./r_ij) ** 12
    pairwise_energy = (r12_term - r6_term)*4
    return pairwise_energy.sum() 

def calculate_pair_energy_np(coordinates, i_particle, box_length, cutoff):
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

    #creates a list of the coordinates for the i_particle
    i_position = coordinates[i_particle]
    
    num_atoms = len(coordinates)
    #i_position has the coordinates for the 1st particle and now we want to see the energy between it and all of its particles
    #that are <cutoff away.
    
    i_position_np = np.array(i_position)
    #new_coords = [x for x in coordinates if x != i_position]
    coordinates_np = np.array(coordinates)
    
    rij = calculate_distance_np(i_position_np, coordinates_np, box_length)
        #should return a list of distances from particle 1 to particle n
        #Now check and see which values are < cutoff and not 0 (distances = 0 are the molecule with itself)
    less_than_cutoff_values = rij[rij<cutoff]
    less_than_cutoff_values3 = less_than_cutoff_values[less_than_cutoff_values != 0]
        #now send all these values to calculate_LJ_np
    e_total = calculate_LJ_np(less_than_cutoff_values3)

    return e_total

def calculate_total_energy_np(coordinates, box_length, cutoff):
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
    #similar to calculate_pair_energy_np however in this version we want the TOTAL energy, not just 1 molecule to all others
    #and so we want all molecules to all others.
    total_energy = 0
    num_atoms = len(coordinates)
    #coordinates_np = np.array(coordinates)
    for i in range(num_atoms):
            
            #The commented lines are unneccesary since we are sending the coordinates of one molecule + all the other molecules
            #thus we can just call the calculate_pair_energy_np function which does the same thing
            
            #mol_position = coordinates[i]
            #mol_position_np = np.array(mol_position)
            
            # Calculate the distance between the particles
            
            #dist_ij = calculate_distance_np(coordinates_np, mol_position_np, box_length)
            #less_than_cutoff_values2 = dist_ij[dist_ij<cutoff]
            #less_than_cutoff_values4 = less_than_cutoff_values2[less_than_cutoff_values2 != 0]
            
            # Add to total energy.
            
            #total_energy += calculate_LJ_np(less_than_cutoff_values4)
            total_energy += calculate_pair_energy_np(coordinates, i, box_length, cutoff)
            
            #now divide total_energy by 2 because we counted the energy for every atom twice
    return total_energy/2
