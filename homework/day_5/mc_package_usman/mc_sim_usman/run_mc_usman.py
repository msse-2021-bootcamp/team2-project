
import math
import random
import numpy as np
import mc_sim_usman.monte_carlo_numpy_usman as mcn
import mc_sim_usman.monte_carlo_usman as mc
import mc_sim_usman.utilities_usman as ut
import mc_sim_usman.time_avg_rdf_usman as rdf

def run_simulation_np(coordinates, box_length, cutoff, reduced_temperature, num_steps, r=5, dr=0.2, max_displacement=0.1, freq=1000, engine=None):
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
    engine: string 
        Use "numpy" to use the numpy implementation, or leave it blank to use the standard python library
    r: int
        Length of the sphere used to calculate the rdf
    dr: float
        Size of the shell of the sphere used to calculate the rdf
    """
    if engine == 'numpy':
        # Calculated quantities
        beta = 1 / reduced_temperature
        num_particles = len(coordinates)
        gr_running_total = np.zeros(r*10-1)

        # Energy calculations
        total_energy = mcn.calculate_total_energy_np(coordinates, box_length, cutoff)
        total_correction = ut.calculate_tail_correction(num_particles, box_length, cutoff)
        total_energy += total_correction


        for step in range(num_steps):
        # 1. Randomly pick one of the particles.
            random_particle = random.randrange(num_particles)
            
            # 2. Calculate the interaction energy of the selected particle with the system.
            current_energy = mcn.calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)
            
            # 3. Generate a random x, y, z displacement.
            x_rand = random.uniform(-max_displacement, max_displacement)
            y_rand = random.uniform(-max_displacement, max_displacement)
            z_rand = random.uniform(-max_displacement, max_displacement)
            
            # 4. Modify the coordinate of Nth particle by generated displacements.
            coordinates[random_particle][0] += x_rand
            coordinates[random_particle][1] += y_rand
            coordinates[random_particle][2] += z_rand
            
            # 5. Calculate the interaction energy of the moved particle with the system and store this value.
            proposed_energy = mcn.calculate_pair_energy_np(coordinates, random_particle, box_length, cutoff)
            delta_energy = proposed_energy - current_energy
            
            # 6. Calculate if we accept the move based on energy difference.
            accept = ut.accept_or_reject(delta_energy, beta)
            
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
                gr_values, r_values = rdf.radial_distribution(coordinates, box_length, r, dr)
                gr_values_np = np.array(gr_values)
                gr_running_total += gr_values_np

    else:
        # Calculated quantities
        beta = 1 / reduced_temperature
        num_particles = len(coordinates)
        gr_running_total = np.zeros(r*10-1)

        # Energy calculations
        total_energy = mc.calculate_total_energy(coordinates, box_length, cutoff)
        total_correction = ut.calculate_tail_correction(num_particles, box_length, cutoff)
        total_energy += total_correction


        for step in range(num_steps):
        # 1. Randomly pick one of the particles.
            random_particle = random.randrange(num_particles)
            
            # 2. Calculate the interaction energy of the selected particle with the system.
            current_energy = mc.calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
            
            # 3. Generate a random x, y, z displacement.
            x_rand = random.uniform(-max_displacement, max_displacement)
            y_rand = random.uniform(-max_displacement, max_displacement)
            z_rand = random.uniform(-max_displacement, max_displacement)
            
            # 4. Modify the coordinate of Nth particle by generated displacements.
            coordinates[random_particle][0] += x_rand
            coordinates[random_particle][1] += y_rand
            coordinates[random_particle][2] += z_rand
            
            # 5. Calculate the interaction energy of the moved particle with the system and store this value.
            proposed_energy = mc.calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
            delta_energy = proposed_energy - current_energy
            
            # 6. Calculate if we accept the move based on energy difference.
            accept = ut.accept_or_reject(delta_energy, beta)
            
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
                gr_values, r_values = rdf.radial_distribution(coordinates, box_length, r, dr)
                gr_values_np = np.array(gr_values)
                gr_running_total += gr_values_np

    time_average_gr = gr_running_total / (num_steps/freq)
     
    return coordinates, time_average_gr, r_values

