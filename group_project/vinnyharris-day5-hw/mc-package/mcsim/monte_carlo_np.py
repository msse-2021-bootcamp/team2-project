from .utils import *
import random
import numpy as np
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


def calculate_lj_np(r_ij):
    """
    The LJ interaction energy between two particles.

    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : np.array of r_ij
        The distance between the particles in reduced units.
    
    Returns
    -------
    pairwise_energy : np.array(float)
        The pairwise Lennard Jones interaction energy in reduced units.

    Examples
    --------
    >>> calculate_LJ(1)
    0

    """
    r1_term = 1 / r_ij
    r6_term = np.power(r1_term, 6)
    r12_term = np.power(r6_term, 2)
    
    pairwise_energy_np = 4 * (r12_term - r6_term)
    
    pairwise_energy = pairwise_energy_np.sum()
    return pairwise_energy


def calculate_total_energy_np(coordinates_np, box_length, cutoff):
    """
    Calculate the total energy of a set of particles using the Lennard Jones potential.
    
    Parameters
    ----------
    coordinates_np : np.array
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
    energy = 0
    
    for i in coordinates_np:
        c_np = np.array(i)
        
        # Get distance
        distance_np = calculate_distance_np(coordinates_np, c_np, box_length)
        
        # Cutoff
        distance = distance_np[distance_np < cutoff]
        
        # Non zero
        distance_nz = distance[distance != 0]
        
        total_energy += calculate_lj_np(distance_nz)

    return total_energy.sum()/2


def calculate_pair_energy_np(coordinates_np, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system)
    
    Parameters
    ----------
    coordinates_np : np.array
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
    i_np = np.array(coordinates_np[i_particle])
    
    distances = calculate_distance_np(i_np, coordinates_np, box_length)
    
    distances_cutoff = distances[distances < cutoff]
    distance_positive = distances_cutoff[distances_cutoff != 0]
    
    e_total = calculate_lj_np(distance_positive)
       
    return e_total


def run_simulation_np(coordinates_np, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000):
    """
    Run a Monte Carlo simulation with the specific parameters and using NumPy for performance gain.

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
    start_time = time.time()
    # Calculated quantities
    beta = 1 / reduced_temperature
    num_particles = len(coordinates_np)

    # Energy calculations
    total_energy = calculate_total_energy_np(coordinates_np, box_length, cutoff)
    #print(total_energy)
    total_correction = calculate_tail_correction(num_particles, box_length, cutoff)
    #print(total_correction)
    total_energy += total_correction


    for step in range(num_steps):
        # 1. Randomly pick one of the particles.
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate the interaction energy of the selected particle with the system.
        current_energy = calculate_pair_energy_np(coordinates_np, random_particle, box_length, cutoff)
        
        # 3. Generate a random x, y, z displacement.
        rand_np=np.random.uniform(-max_displacement, max_displacement, 3)
        
        # 4. Modify the coordinate of Nth particle by generated displacements.
        coordinates_np[random_particle]= coordinates_np[random_particle] + rand_np
        
        # 5. Calculate the interaction energy of the moved particle with the system and store this value.
        proposed_energy = calculate_pair_energy_np(coordinates_np, random_particle, box_length, cutoff)
        delta_energy = proposed_energy - current_energy
        
        # 6. Calculate if we accept the move based on energy difference.
        accept = accept_or_reject(delta_energy, beta)

        # 7. If accepted, move the particle.
        if accept:
            total_energy += delta_energy
        else:
            #Move not accepted, roll back coordinates
            coordinates_np[random_particle]= coordinates_np[random_particle] - rand_np
        
        # 8. Print the energy if step is a multiple of freq.
        #if step % freq == 0:
        #    print(step, total_energy/num_particles)
    
    end_time = time.time()
    elapsed_time = (end_time-start_time)
    return elapsed_time