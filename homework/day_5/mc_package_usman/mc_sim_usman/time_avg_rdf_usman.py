import math
import random
import mc_sim_usman as mcn

def radial_distribution(coordinates, box_length, r, dr):
    """
    The radical distribution analysis of a box containing n particles after n steps of a Monte Carlo simulation.
    
    Calculates the pair correlation function g(r) at as many r as you want for a fixed dr.
    
    Parameters
    ----------
    coordinates : list
        The final coordinates for all the particles in the system after n MC steps.
            
    box_length : float
        The length of the box for periodic bounds used in the calculate_distance function.
    
    r : int
        The length of the distance from one particle (size of the shell).
    
    dr : float
        The length of the distance interval (size of the shell thickness).
        
    Returns
    -------
    gr : list of float
        The list of pair correlation function g(r) for each r we investigated.
        
    rlist : list of float
        The list of r/density for each r which we investigated.
    """
    gr=[]
    rlist=[]
    volume = box_length**3
    num_particles = len(coordinates)
    density = num_particles/volume
    for i in range(1, r*10):
        i = i*0.1
        count = 0
        for particle_1 in range(num_particles):
            for particle_2 in range(num_particles):
                if particle_1 != particle_2:
                    dist_between_particles = mcn.calculate_distance(coordinates[particle_1], coordinates[particle_2], box_length)
                    if dist_between_particles > i and dist_between_particles < i+dr:
                        count += 1
                        
        div_count = count/num_particles
        div_sphere_volume = div_count/(4*math.pi*i**2*dr)
        div_particle_num_density = div_sphere_volume/(density)
        rlist.append(i/density)
        gr.append(div_particle_num_density)
    return gr, rlist