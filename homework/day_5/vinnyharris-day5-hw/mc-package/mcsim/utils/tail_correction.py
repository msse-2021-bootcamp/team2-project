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