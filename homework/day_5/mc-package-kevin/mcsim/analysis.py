"""
Pressure analysis of MC 
"""

def calculate_P_virial(r_ij):
    """
    The virial between two particles.

    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.
    
    Returns
    -------
    p_virial : float
        The pairwise Lennard Jones interaction energy in reduced units.

    Examples
    --------
    >>> calculate_LJ(1)
    0

    """
    
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    
    p_virial = -48 * (r12_term - r6_term)
    
    return p_virial

def calculate_pressure(p_virial,num_particles, temperature,box_length):

    volume = box_length**3
    pressure = (3*num_particles*temperature- p_virial)/(3*volume)
    
    return pressure

def calculate_pressure_tail_correction(num_particles, box_length, cutoff):
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
    brackets = (2/3*math.pow(1/cutoff,9)) - math.pow(1/cutoff,3)
    volume = box_length**3
    
    constant = (16*math.pi*(num_particles**2))/(3*math.pow(volume,2))
    
    tail_correction = constant * brackets
    
    return tail_correction

def calculate_P_virial_np(r_ij):
    """
    The virial between two particles.

    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : np.array
        The distance between the particles in reduced units.
    
    Returns
    -------
    p_virial : float
        The pairwise Lennard Jones interaction energy in reduced units.

    """
    r6_term = (1./r_ij)** 6
    r12_term = (1./r_ij)** 12
    
    p_virial = -48 * (r12_term - r6_term)
    
    return p_virial




