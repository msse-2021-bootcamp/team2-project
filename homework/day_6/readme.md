Homework for day 6

This program evaluates the integrand 1 / (1 + X^2) with the bounds of a,b on n_points

1. The function to evaluate the integrand will receive a parameter of type double and will return a double as well.

2. The function to calculate the approximated value of the integral will include a parameter list with the bounds of the integral a,b as type double, also as part of the param list we will include the nPoints type int and for the return type will be double.

3. The width of the rectangles will evaluate to (b-a) / nPoints

4. For every increment of dx, we need to find the midpoint between the current and the next dx. 
Where dx is defined as (b-a) / nPoints, for all points i (nPoints) from a to b every idx + dx/2
The other way would be as (idx + (1+i)dx ) / 2

5. The integration function calculates the total area by summing all the rectangle areas. 
For width we are using dx = (b-a) / nPoints
For height fx = 1 / (1 + X^2)

rectangleArea = fx * dx

Sum of all rectangleAreas

6. The main function handles the call to the integration function and printing the result.
Based on our calculations we observed that 40 nPoints give a good 4 decimal approximation.



