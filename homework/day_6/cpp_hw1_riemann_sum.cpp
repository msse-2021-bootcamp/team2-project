#include <iostream>

double eval_integrand(double x)
{
    return (1/(1+(x*x)));
}

double integration_function(int a, int b, int npoints)
{
    double integeral_sum=0;
    float domain = ((float)b-a);
    float width = ((float)b-a)/npoints;

    for (float i=a; i < b; i += width){
        integeral_sum += eval_integrand(i+(width/2));
        double midpoint = i+(width/2);
        std::cout << midpoint<< std::endl;
    }
    return ((domain/(npoints))*integeral_sum);
}

int main()
{
    std::cout <<"The function integrates to: "<< integration_function(0, 1, 10000) << " at 10,000 steps!" << std::endl;
    return 0;
}

// Answers to questions:
// The eval_integrand function should receive a double and return a double
// a, b and npoints should all be integers
// dx=(b-a)/npoints is the width of the rectangles we will be using
// we will need to evaluate the integrand at a+dx/2 for each rectangle and move a by +dx every increment
//after summing all the y values we can multiply the integral sum by (domain/npoints) which is the same as calculating the area of each rectangle during the for loop
//and summing them up.
