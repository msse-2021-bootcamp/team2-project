"""
Functions for running a Motne Carlo Simulation using numpy arrays
"""

import math
import random
import numpy as np

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

def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.
    
    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.
       
    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates
    """
    
    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()
    
    atomic_coordinates = []
    
    for atom in coordinates:
        split_atoms = atom.split()
        
        float_coords = []
        
        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))
            
        atomic_coordinates.append(float_coords)
    atomic_coordinates_np=np.array(atomic_coordinates)       
    return atomic_coordinates_np, box_length

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

def calculate_tail_correction(num_particles, box_length, cutoff):
    """
    The tail correction associated with using a cutoff radius.
    
    Computes the tail correction based on a cutoff radius used in the LJ energy calculation in reduced units.
    
    Parameters
    ----------
    num_particles : int
        The number of particles in the system.
    
    box_length : int
        Size of the box length of the system, used to calculate volume.
    
    cutoff : int
        Cutoff distance.
    
    Returns
    -------
    tail_correction : float
        The tail correction associated with using the cutoff.
    """
    
    brackets = (1/3*math.pow(1/cutoff,9)) - math.pow(1/cutoff,3)
    volume = box_length**3
    
    constant = ((8*math.pi*(num_particles**2))/(3*volume))
    
    tail_correction = constant * brackets
    
    return tail_correction

def accept_or_reject(delta_U, beta):
    """
    Accept or reject a move based on the Metropolis criterion.
    
    Parameters
    ----------
    detlta_U : float
        The change in energy for moving system from state m to n.
    beta : float
        1/temperature
    
    Returns
    -------
    boolean
        Whether the move is accepted.
    """
    if delta_U <= 0.0:
        accept = True
    else:
        #Generate a random number on (0,1)
        random_number = random.random()
        p_acc = math.exp(-beta*delta_U)
        
        if random_number < p_acc:
            accept = True
        else:
            accept = False
    return accept

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

def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement=0.1, freq=1000):
    """
    Run a Monte Carlo simulation with the specific parameters.

    Parameters
    ---------
    coordinates: list
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
    total_energy = calculate_total_energy(coordinates, box_length, cutoff)
    total_correction = calculate_tail_correction(num_particles, box_length, cutoff)
    total_energy += total_correction


    for step in range(num_steps):
    # 1. Randomly pick one of the particles.
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate the interaction energy of the selected particle with the system.
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        # 3. Generate a random x, y, z displacement.
        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)
        
        # 4. Modify the coordinate of Nth particle by generated displacements.
        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand
        
        # 5. Calculate the interaction energy of the moved particle with the system and store this value.
        proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        delta_energy = proposed_energy - current_energy
        
        # 6. Calculate if we accept the move based on energy difference.
        accept = accept_or_reject(delta_energy, beta)
        
        # 7. If accepted, move the particle.
        if accept:
            total_energy += delta_energy
        else:
            #Move not accepted, roll back coordinates
            coordinates[random_particle][0] -= x_rand
            coordinates[random_particle][1] -= y_rand
            coordinates[random_particle][2] -= z_rand
        
        # 8. Print the energy if step is a multiple of freq.
        if step % freq == 0:
            print(step, total_energy/num_particles)

    return coordinates

print("package imported!")