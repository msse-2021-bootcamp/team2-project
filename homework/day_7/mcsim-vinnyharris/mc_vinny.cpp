#include <iostream>
#include <random>       // for random numbers
#include <chrono>       // for generating the seed
#include <fstream>      // for reading/writing files
#include <array>        // for std::array
#include <vector>       // for std::vector
#include <utility>      // for std::pair
#include <cmath>        // for std:math
#include <time.h>       // for time

using std::cout;
using std::endl;
using std::string;
using std::pow;

// Make some types more convenient
typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

/*! Generate a random double within a given range */
double randomDouble(double lower_bound, double upper_bound)
{
   std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
   return dist(re);
}

/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double randomInteger(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
}

std::pair<Coordinates, double> read_xyz(string file_path)
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


/*
    Calculates the distance between two points. When box_length is set, 
    the minimum image convention is used to calculate the distance between the points.

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
*/

double calculateDistance(AtomCoord coord1, AtomCoord coord2, double boxLength)
{
    double distance = 0.0;
    int count = coord1.size();

    for (size_t i = 0; i < count; i++)
    {
        double holdDist = abs(coord2.at(i) - coord1.at(i));

        if(boxLength) 
        {
            holdDist -= (boxLength * round(holdDist/boxLength));
        }

        distance += (holdDist * holdDist);

    }

    return sqrt(distance);
}

/*
    The LJ interaction energy between two particles.
    Computes the pairwise Lennard Jones interaction energy based on the separation 
    distance in reduced units.

    Parameters
    ----------
    rIJ : double
        The distance between the particles in reduced units.
    
    Returns
    -------
    pairwise_energy : double
        The pairwise Lennard Jones interaction energy in reduced units.
*/
double calculateLJ(double rIJ)
{   
    double r6 = pow(1/rIJ, 6.0);
    double r12 = pow(r6, 2.0);

    return 4 * (r12 - r6);
}

/*
    Calculate the total energy of a set of particles using the Lennard Jones potential.
    
    Parameters
    ----------
    coords : list
        A nested list containing the x, y,z coordinate for each particle
    boxLength : double
        The length of the box. Assumes cubic box.
    cutoff : double
        The cutoff length
    
    Returns
    -------
    total_energy : double
        The total energy of the set of coordinates.
*/
double calculateTotalEnergy(Coordinates coords, double boxLen, double cutoff)
{
    double totalEnergy = 0.0;
    int nAtoms = coords.size();
    cout << "nAtoms: " << nAtoms << endl;

    for (size_t i = 0; i < (nAtoms-1); i++)
    {
        for (size_t j = i+1; j < nAtoms; j++)
        {
            // # Calculate the distance between the particles
            double distIJ = calculateDistance(coords.at(i), coords.at(j), boxLen);

            if(distIJ < cutoff) {
                double ljIj = calculateLJ(distIJ); 
                
                // Add to total energy.
                totalEnergy += ljIj;
            }
            
        } 
    }
    return totalEnergy;
}

/*
    Calculates the interaction energy of a particle with its environment 
    (all other particles in the system).
    
    Parameters
    ----------
    coords : vector
        The coordinates for all the particles in the system.
        
    iParticle : int
        The particle number for which to calculate the energy. Index 0 to n.
        
    cutoff : double
        The simulation cutoff. Beyond this distance, interactions are not calculated.
    
    boxLength : double
        The length of the box for periodic bounds
        
    Returns
    -------
    eTotal : double
        The pairwise interaction energy of the ith particles with all other particles 
        in the system.
    
*/
double calculatePairEnergy(Coordinates coords, int iParticle, double boxLength, double cutoff) 
{
    double eTotal = 0.0;
    AtomCoord iAtom = coords[iParticle];
    int numAtoms = coords.size();

    for (size_t i = 0; i < numAtoms; i++)
    {
        if(i == iParticle) { continue; }
        
        AtomCoord currentAtom = coords[i];

        double rij = calculateDistance(iAtom, currentAtom, boxLength);
        if(rij < cutoff) {
            double lj = calculateLJ(rij);
            eTotal += lj;
        }
    }

    return eTotal;
}


/*
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
*/

bool acceptOrReject(double deltaU, double beta)
{   
    bool accept = false;
    if(deltaU <= 0.0) {
        accept = true;
    }    
    else {
        //Generate a random number on (0,1)
        double randomNum = randomDouble(0.0, 1.0);
        double pAcc = exp((-beta)*deltaU);
        
        if(randomNum < pAcc) {
            accept = true;
        }    
        else { 
            accept = false;
        }    
    } 
    return accept;
}

/*
    Calculates the tail correction
    
    Parameters
    ----------
    box_legth : double
        The distance between the particles in reduced units.
    n_particles : int
        The number of particles.
    cutoff : double
        The cutoff value.
    
*/
double calculateTailCorrection(double boxLen, int nParticles, double cutoff)
{
    double piTerm = (8.0 * M_PI * pow(nParticles, 2.0)) / (3.0 * pow(boxLen, 3.0));
    double r9Term = (1.0/3.0) * (pow(1.0/cutoff, 9.0));
    double r3Term = pow(1.0/cutoff, 3.0);

    return piTerm * (r9Term - r3Term);
}

/*
    Runs Monte Carlo Simulation
*/
void runSimulation(Coordinates coords, double boxLen, double cutoff, double reducedTemp, int numSteps, double maxDisplacement, int frequency=1000)
{   
    // Calculated quantities
    double beta = 1.0 / reducedTemp;
    int numParticles = coords.size();

    // Energy calculations
    double totalEnergy = calculateTotalEnergy(coords, boxLen, cutoff);
    double totalCorrection = calculateTailCorrection(boxLen, numParticles, cutoff);

    totalEnergy += totalCorrection;

    for (size_t i = 0; i < numSteps; i++)
    {
        // 1. Randomly pick one of the particles.
        double randomParticle = randomInteger(0, numParticles);

        // 2. Calculate the interaction energy of the selected particle with the system.
        double currentEnergy = calculatePairEnergy(coords, randomParticle, boxLen, cutoff);
        
        // 3. Generate a random x, y, z displacement.
        double xRand = randomDouble(-maxDisplacement, maxDisplacement);
        double yRand = randomDouble(-maxDisplacement, maxDisplacement);
        double zRand = randomDouble(-maxDisplacement, maxDisplacement);

        // 4. Modify the coordinate of Nth particle by generated displacements.
        coords[randomParticle][0] += xRand;
        coords[randomParticle][1] += yRand;
        coords[randomParticle][2] += zRand;

        // 5. Calculate the interaction energy of the moved particle with the system and store this value.
        double proposedEnergy = calculatePairEnergy(coords, randomParticle, boxLen, cutoff);
        double deltaEnergy = proposedEnergy - currentEnergy;
        
        // 6. Calculate if we accept the move based on energy difference.
        bool accept = acceptOrReject(deltaEnergy, beta);

        // # 7. If accepted, move the particle.
        if (accept) {
            totalEnergy += deltaEnergy;
        } 
        else {
            // move not accepted, roll back coordinates
            coords[randomParticle][0] -= xRand;
            coords[randomParticle][1] -= yRand;
            coords[randomParticle][2] -= zRand;
        }

        // 8. Print the energy if step is a multiple of freq.
        if(i % frequency == 0) {
            //print(step, total_energy/num_particles)
            cout << i << " " <<  totalEnergy/numParticles << endl;
        }
                   
    }
    
}

/*
    Validates functions in this file
*/
void validate() {

    double cutoff = 3.0;
    double boxLen = 10.0;
    AtomCoord c1 = {0, 0, 0};
    AtomCoord c2 = {0, 0, 8};

    cout << "Distance: " << calculateDistance(c1, c2, boxLen) << ", expected is 2" << endl;

    double rIJ = pow(2.0, (1.0/6.0));
    cout << "LJ: " << calculateLJ(rIJ) << ", expected: -1" << endl;

    AtomCoord a1 = {0.0, 0.0, 0.0};
    AtomCoord a2 = {0.0, 0.0, rIJ};
    AtomCoord a3 = {0.0, 0.0, 2.0 * rIJ};
    
    Coordinates coords;
    coords.push_back(a1);
    coords.push_back(a2);
    coords.push_back(a3);
    
    cout << "Total energy: " << calculateTotalEnergy(coords, boxLen, 3.0) << ", expected: -2" << endl;

    double tailC = calculateTailCorrection(boxLen, 800, cutoff);
    double tailExpected = -198.4888837441566;

    cout << "Tail correction: " << tailC << ", expected: " << tailExpected << endl;

    cout << "Accept or reject: " << acceptOrReject(1.0, 1.0) << ", expected: true" << endl;

}


/*
    Main
*/
int main(void)
{   
    // Initialize random number generation based on time
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
    
    validate();

    std::pair<Coordinates, double> xyzInfo = read_xyz("../../../../lj_sample_configurations/lj_sample_config_periodic1.txt");

    Coordinates coordinates= xyzInfo.first;
    double boxLength = xyzInfo.second;

    double cutoff = 3.0;
    int numSteps = 50000;
    double reducedTemperature = 0.9;
    double maxDisplacement = 0.1;
    int freq = 1000;

    // Start measuring time
    time_t start, end;
    time(&start);
    runSimulation(coordinates, boxLength, cutoff, reducedTemperature, numSteps, maxDisplacement, freq);
    time(&end);
    time_t elapsedTime = end - start;
 
    cout << "Elapsed time: " << elapsedTime << " secomds" << endl;

    return 0;
}