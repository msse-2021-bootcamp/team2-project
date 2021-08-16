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