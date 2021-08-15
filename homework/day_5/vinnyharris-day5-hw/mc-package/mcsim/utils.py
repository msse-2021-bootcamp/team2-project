import math
import random


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

def calculate_tail_correction(box_legth, n_particles, cutoff):
    """
    Calculates the tail correction
    
    Parameters
    ----------
    box_legth : float
        The distance between the particles in reduced units.
    n_particles : int
        The number of particles.
    cutoff : float
        The cutoff value.
    """
    pi_term = (8 * math.pi * math.pow(n_particles, 2)) / (3 * math.pow(box_legth, 3))
    r9_term = (1/3) * (math.pow(1/cutoff, 9))
    r3_term = math.pow(1/cutoff, 3)
    
    tail_correction = pi_term * (r9_term - r3_term)
    
    return tail_correction


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
        
    return atomic_coordinates, box_length