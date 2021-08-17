#include<iostream>

double integrand(double x)
{
    return 1/(1+x*x);
}

double integral(double a, double b, int n_points)
{
    double Riemann = 0 ;
    double dx = (b - a)/n_points;
    for (int i=0; i<n_points; i++)
        {
            double boxArea = dx*integrand(a + i*dx + dx/2);
            Riemann += boxArea;
        }
    return Riemann;
}

int main (void)
{
    for (int i=39; i<190;i++)
        {
            double pi = 4*integral(0,1,i);
            std::cout<<"For "<<i+1<<" boxes, pi is approximately " << pi << "."<<std::endl;
        }
    return 0;
}