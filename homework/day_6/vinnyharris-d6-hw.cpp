

#include <iostream>

double calculateFx(double x) {

    return 1 / (1 + x * x);
}


double calculateIntegral(double a, double b, int nPoints) {

    // Calculate dx
    double dx = (b-a) / nPoints;
    double totalArea = 0.0;

    for(int i = 0; i < nPoints; i++) {

        // get checkPoint
        double midPoint = i * dx + dx/2.0;

        //calculate integrant
        double fx = calculateFx(midPoint);

        // Get rectangle area
        totalArea += dx * fx;

    }

    return totalArea;
}

int main() {
    double a = 0.0;
    double b = 1.0;
    int nPoints = 40;

    double approximation = calculateIntegral(a, b, nPoints);
    std::cout << "Approximation: " << approximation << "\n";
    std::cout << "Pi: " << (approximation * 4);
    
}

