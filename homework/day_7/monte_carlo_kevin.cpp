#include <iostream>
#include <fstream> // for reading/writing files
#include <array>   // for std::array
#include <vector>  // for std::vector
#include <utility> // for std::pair
#include <random> // for random numbers
#include <chrono> // for generating the seed
#define _USE_MATH_DEFINES
#include <cmath>   // for std::pow, std::sqrt, std::rint

typedef std::array <double,3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

double calculateDistance (AtomCoord coord1, AtomCoord coord2, double boxLength)
{   
    double dSquared = 0;
    for(int i=0;i<3;i++)
    {   //Calculate distance 
        double dimDist =coord2[i]-coord1[i];
        if (boxLength>0) //Apply periodic boundary condition if given
            dimDist = dimDist - boxLength*std::rint(dimDist/boxLength);
        dSquared = dSquared + std::pow(dimDist,2);
    }
    double distance = std::sqrt(dSquared);
    return distance;
}

double calculateLJ(double rij)
{
    double r6Term = std::pow(1.0/rij,6.0);
    double r12Term = std::pow(r6Term,2.0);
    double pairwiseEnergy = 4.0* (r12Term-r6Term);
    return pairwiseEnergy;
}

double calculatePairEnergy(Coordinates coordinates, int particleIndex, double boxLength, double cutoff)
{
    double pairEnergy=0.0;
    AtomCoord iCoords = coordinates[particleIndex];
    int numParticles = coordinates.size();
    for (int j=0;j<numParticles;j++)
    {
    if (particleIndex != j)
        {
        AtomCoord jCoords = coordinates[j];
        double rij = calculateDistance(iCoords,jCoords,boxLength);
        if (rij<cutoff)
        {
            double ePair = calculateLJ(rij);
            pairEnergy += ePair;
        }}
    }
    return pairEnergy;
}

double calculateTotalEnergy(Coordinates coordinates, double boxLength, double cutoff)
{
    double totalEnergy = 0;
    int numParticles = coordinates.size();
    for (int index=0;index<numParticles;index++)
    {
        double pairEnergy=calculatePairEnergy(coordinates, index, boxLength, cutoff);
        totalEnergy += pairEnergy/2.0;
    }
    return totalEnergy;
}

double CalculateTailCorrection(int numParticles, double boxLength, int cutoff)
{
    double constant = (8.0*M_PI*(std::pow(numParticles,2.0)))/(3.0*std::pow(boxLength,3.0));
    double brackets = (1.0/3.0*std::pow(1.0/cutoff,9.0)) - std::pow(1.0/cutoff,3.0);
    double tailCorrection = constant * brackets;
    return tailCorrection;
}

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

// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

/*! Generate a random double within a given range */
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

bool acceptORreject(double deltaE,double beta)
{
    bool accept;
    if (deltaE<0)
        {accept=true;}
    else
    {
        double random=random_double(0,1);
        double pAcc = std::exp(-beta*deltaE);
        if (random<pAcc)
            {accept=true;}
        else
            {accept=false;}
    return accept;
    }
}

void runMCcpp (Coordinates coordinates, double cutoff, double reducedTemp, int numSteps, double max_displacement=0.1, int freq=1000, double boxLength = 0)
{
    // Calculated quantities
    double beta = 1.0/ reducedTemp;
    int numParticles = coordinates.size();
    // Energy calculations
    double totalEnergy = calculateTotalEnergy(coordinates,boxLength,cutoff);
    double tailCorrection = CalculateTailCorrection(numParticles,boxLength,cutoff);
    totalEnergy += tailCorrection;
    for (int step=0;step<numSteps;step++)
    {
        int randomParticle = random_integer(0,numParticles);
        double currentEnergy = calculatePairEnergy(coordinates,randomParticle,boxLength,cutoff);
        
        double randX = random_double(-max_displacement,max_displacement);
        double randY = random_double(-max_displacement,max_displacement);
        double randZ = random_double(-max_displacement,max_displacement);

        coordinates[randomParticle][0] += randX;
        coordinates[randomParticle][1] += randY;
        coordinates[randomParticle][2] += randZ;
        
        double proposedEnergy = calculatePairEnergy(coordinates,randomParticle,boxLength,cutoff);
        double deltaEnergy = proposedEnergy - currentEnergy;
        bool accept = acceptORreject(deltaEnergy,beta);
        
        if (accept)
            {totalEnergy += deltaEnergy;}
        else
        {
            coordinates[randomParticle][0] -= randX;
            coordinates[randomParticle][1] -= randY;
            coordinates[randomParticle][2] -= randZ;
        }
        if (step % freq == 0)
        {
            std::cout<<"At step "<<step<<", the total Energy per particle is "<<totalEnergy/numParticles<<std::endl;
        }
    }
}

int main(void)
{
    std::pair<Coordinates, double> xyz_info = read_xyz("../../lj_sample_configurations/lj_sample_config_periodic1.txt");

    Coordinates coords = xyz_info.first;
    double boxLength = xyz_info.second;
    double cutoff = 3;
    double reduceTemp = 0.9;
    int steps = 50000;
    int freq = 1000;
    // Initialize random number generation based on time
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    runMCcpp(coords, cutoff, reduceTemp, steps, boxLength);

    return 0;
}