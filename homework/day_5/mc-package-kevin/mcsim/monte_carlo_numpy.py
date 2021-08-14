"""
Functions for running a Monte Carlo Simulation using numpy arrays
"""

import math
import random
import numpy as np
import mcsim.utilities as util

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
    distance : float
        The distance between the two points accounting for periodic boundaries
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

#calculate LJ for num py
def calculate_LJ_np(r_ij):
    """
    The LJ interaction energies for an array of distances.
    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units. also considers cutoff

    Parameters
    ----------
    r_ij : np.arry
        The distance between the particles in reduced units.

    Returns
    -------
    pairwise_energy : np.array
        The pairwise Lennard Jones interaction energy in reduced units.
    """
    r6_term = (1./r_ij)** 6
    r12_term = (1./r_ij)** 12
    pairwise_energy = 4 * (r12_term - r6_term)

    return pairwise_energy

def calculate_pair_energy_np(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system)
    
    Parameters
    ----------
    coordinates : np.array
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
    
    i_position_np =np.array(coordinates[i_particle])
    distances= calculate_distance_np(i_position_np,coordinates,box_length)
    #only account for distances less than cut off
    distances=distances[distances<cutoff]
    #eliminates point interaction with self
    new_distances=distances[distances != 0]
    e_pair = calculate_LJ_np(new_distances)
    e_total = e_pair.sum()
    
    return e_total

def calculate_total_energy_np(coordinates, box_length, cutoff):
    """
    Calculate the total energy of a set of particles using the Lennard Jones potential.
    
    Parameters
    ----------
    coordinates : np.array
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
    pair_energies=[]
    for i in range(coordinates.shape[0]):
        #calculate_pair_energy_np will eliminate for cutoff and self interaction
        particle_pair_energy=calculate_pair_energy_np(coordinates, i, box_length, cutoff)
        pair_energies.append(particle_pair_energy)
    #turn list of pair energies into a numpy array
    pair_energies_np=np.array(pair_energies)
    #sum all pair energies to calculate total energy
    sum_energy=pair_energies_np.sum()
    # all pair-wise interactions are accounted for twice in the model
    total_energy=sum_energy/2
    return total_energy


def run_simulation_np(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000):
    """
    Run a Monte Carlo simulation with the specific parameters implementing numpy arrays.

    Parameters
    ---------
    coordinates: np.array
    box_length: float
    cutoff: float
    reduced_temperature: float
    num_steps:int
    max_displacement: float
    freq_checks: int
    
    """
    # Calculated quantities
    beta = 1 / reduced_temperature
    num_particles = len(coordinates)

    # Energy calculations
    total_energy = calculate_total_energy_np(coordinates, box_length, cutoff)
    total_correction = util.calculate_tail_correction(num_particles, box_length, cutoff)
    total_energy += total_correction


    for step in range(num_steps):
        # 1. Randomly pick one of the particles.
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate the interaction energy of the selected particle with the system.
        current_energy = calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)
        
        # 3. Generate a random x, y, z displacement.
        rand_np=np.random.uniform(-max_displacement, max_displacement, 3)
        
        # 4. Modify the coordinate of Nth particle by generated displacements.
        coordinates[random_particle]= coordinates[random_particle] + rand_np
        
        # 5. Calculate the interaction energy of the moved particle with the system and store this value.
        proposed_energy = calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)
        delta_energy = proposed_energy - current_energy
        
        # 6. Calculate if we accept the move based on energy difference.
        accept = util.accept_or_reject(delta_energy, beta)
        
        # 7. If accepted, move the particle.
        if accept:
            total_energy += delta_energy
        else:
            #Move not accepted, roll back coordinates
            coordinates[random_particle]= coordinates[random_particle] - rand_np
        
        # 8. Print the energy if step is a multiple of freq.
        if step % freq == 0:
            print(step, total_energy/num_particles)

    return coordinates

print("package imported!")