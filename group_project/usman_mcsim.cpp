#include <fstream> // for reading/writing files
#include <utility> // for std::pair
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random> // for random numbers
#include <chrono> // for generating the seed
#include <fstream> 

using namespace std;
typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;
std::default_random_engine re;


//reads information from the file.
std::pair<Coordinates, double> read_xyz(std::string file_path)
{
    // Opens up a file stream for input
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if(!infile.is_open())
    {   
        throw std::runtime_error("File path in read_xyz does not exist!");
    }
    
    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;
    
    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;
    
    // now the number of atoms
    infile >> num_atoms;
    
    // Uncomment to help troubleshoot
    //std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;
    
    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;
    
    for(int i = 0; i < num_atoms; i++)
    {   
        AtomCoord coord;
        
        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];
        
        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}

double random_double(double lower_bound, double upper_bound)
{
   std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
   return dist(re);
}

/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double random_integer(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
}  

double calculate_lj(double r_ij)
{
    double division = 1.0/r_ij;
    double r6_term = pow(division,  6.0);
    double r12_term = pow(division,  12.0);

    double pairwise_energy = 4.0*(r12_term - r6_term);
    return pairwise_energy;
}

double calculate_distance(AtomCoord mol1, AtomCoord mol2, double box_length)
{
    double distance = 0.0;
    for (int i = 0; i<3; i++)
    {
        double hold_dist = mol1[i] - mol2[i];
        if (box_length > 0.0)
        {
            hold_dist = hold_dist - (box_length * round(hold_dist/box_length));
        } 
        distance += pow(hold_dist, 2.0);
    }
    return sqrt(distance);
}

double calculate_pair_energy(Coordinates coor_all, int mol_position, double box_length, int cutoff)
{
    double e_total = 0.0;
    AtomCoord mol = coor_all[mol_position];
    int num_atoms = coor_all.size();

    for(int i=0; i < num_atoms; i++)
    {
        if (i != mol_position)
        {
            AtomCoord mol2 = coor_all[i];
            double rij = calculate_distance(mol, mol2, box_length);
            if (rij < cutoff)
            {
                double e_pair = calculate_lj(rij);
                e_total += e_pair;
            }
        }
    }
    return e_total;
}

double calculate_total_energy(Coordinates coor_all, double box_length, int cutoff)
{
    int num_atoms = coor_all.size();
    double calculate_total_pair_energy = 0.0;
    for(int i=0; i < num_atoms; i++)
    {
        calculate_total_pair_energy += calculate_pair_energy(coor_all, i, box_length, cutoff);
    }
    return calculate_total_pair_energy/2;
}

double calculate_tail_correction(int num_particles, double box_length, int cutoff)
{
    double brackets = (1.0/3.0*pow(1.0/cutoff,9.0))-pow(1.0/cutoff, 3.0);
    double volume = pow(box_length,3.0);
    double constant ((8.0*M_PI*(pow(num_particles,2.0)))/(3.0*volume));
    double tail_correction = constant * brackets;
    return tail_correction;
}

bool accept_or_reject(double delta_U, double beta)
{
    bool accept;
    if (delta_U <= 0.0)
    {
        accept=true;
    }
    else
    {
        double random_number = random_double(0.0, 1.0);
        double p_acc = exp(-beta*delta_U);
        if (random_number < p_acc)
        {
            accept=true;
        }
        else
        {
            accept=false;
        }
    }
    return accept;
}

void mc_simulation(Coordinates coords, double box_length, int cutoff, double reduced_temperature, int num_steps, double max_displacement, int freq, ofstream & out_file)
{
    auto start = chrono::high_resolution_clock::now();
    int num_particles = coords.size();
    double beta = 1/reduced_temperature;

    //Total energy calculations
    double total_energy = calculate_total_energy(coords, box_length, cutoff);
    //cout << "Initial Total Energy: " << total_energy << endl;
    double total_correction = calculate_tail_correction(num_particles, box_length, cutoff);
    //cout << "Total Correction: " << total_correction << endl;

    total_energy += total_correction;

    for(int i=0; i<num_steps; i++)
    {
        //randomly picks a particle
        int random_particle = random_integer(0, num_particles);

        //calculates the interaction energy of the particle with the system
        double current_energy = calculate_pair_energy(coords, random_particle, box_length, cutoff);

        //generate random x,y,z displacement
        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        //modify the coordinates of the Nth particle by generated displacements
        coords[random_particle][0] += x_rand;
        coords[random_particle][1] += y_rand;
        coords[random_particle][2] += z_rand;

        //calculate the interaction energy of the moved particle with the system and store the value
        double proposed_energy = calculate_pair_energy(coords, random_particle, box_length, cutoff);
        double delta_energy = proposed_energy-current_energy;

        bool accept = accept_or_reject(delta_energy, beta);

        //if move is accepted, move the particle
        if (accept)
        {
            total_energy += delta_energy;
        }
        //move not accepted, roll the coordinates back
        else
        {
        coords[random_particle][0] -= x_rand;
        coords[random_particle][1] -= y_rand;
        coords[random_particle][2] -= z_rand;
        }
        
        /*if (i % freq == 0)
        {
        cout << "Step: " << i << " Average energy: " << total_energy/num_particles << endl;
        }*/
    
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    out_file << " elapsed time  " << elapsed.count() << " s"<<endl;
}
 
int main(void)
{   
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::pair<Coordinates, double> xyz_info = read_xyz("../lj_sample_configurations/lj_sample_config_periodic1.txt");

    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;
    ofstream out_file("MCsim_cpp.txt");

    double cutoff=3.0;
    double reduced_temp = 0.9;
    int steps = 50000;
    double max_displacement = 0,1

    out_file<<"For base case:";
    mc_simulation(coords, box_length, cutoff, reduced_temp, steps, max_displacement, 10000, out_file);
    double temperatures[3] = {0.4, 1.4, 3.1};
    out_file<<"Testing Temperature:"<<endl;
    for (int i=0; i<3; i++)
    {
        for (int trial=1; trial<4; trial++)
        {
            out_file<<"For reduced_temp "<<temperatures[i]<<" , trial "<<trial;
            mc_simulation(coords, box_length, cutoff, temperatures[i], steps, max_displacement, 10000, out_file);}
        }
    double cutoffs[3] = {2.5, 3.5, 4.0};
    out_file<<"Testing cutoff:"<<endl;
    for (int i=0; i<3; i++)
    {
        for (int trial=1; trial<4; trial++)
        {
            out_file<<"For cutoff "<<cutoffs[i]<<" , trial "<<trial;
            mc_simulation(coords, box_length, cutoffs[i], reduced_temp, steps, max_displacement,10000, out_file);}
        }    
    int nums_steps[3] = {10000,30000,70000};
    out_file<<"Testing Number of Steps:"<<endl;
    for (int i=0; i<3; i++)
    {
        for (int trial=1; trial<4; trial++)
        {
            out_file<<"For "<<nums_steps[i]<<" steps, trial "<<trial;
            mc_simulation(coords, box_length, cutoff, reduced_temp, nums_steps[i], max_displacement, 10000, out_file);}
        } 
    double displacements[3] = {0.05,0.15,0.2};
    out_file<<"Testing Max Displacements:"<<endl;
    for (int i=0; i<3; i++)
    {
        for (int trial=1; trial<2; trial++)
        {
            out_file<<"For max displacement of "<<displacements[i]<<" , trial "<<trial;
            mc_simulation(coords, box_length, cutoff, reduced_temp, steps, displacements[i], 10000, out_file,);
        } 
    }
    return 0;

}