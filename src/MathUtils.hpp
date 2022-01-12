#ifndef CAJETE_MATH_UTILS_HPP
#define CAJETE_MATH_UTILS_HPP 

#include <iostream>
#include <math.h>

namespace Cajete 
{

template <typename PositionType, std::size_t N>
double calculate_distance(PositionType (&p1)[N], PositionType (&p2)[N])
{
    double l = 0.0;

    for(auto i = 0; i < N; i++) {
        double diff = p1[i] - p2[i];
        l += diff*diff;
    }
    return sqrt(l);
}

double sigmoid(double input, double coeffcient)
{
    return 1.0/(1.0+exp(-coeffcient*input));
}

double heaviside(double input, double threshold)
{
    if(input >= threshold)
    {
        return 1.0;
    } 
    else
    {
        return 0.0;
    }
}

} //end namespace Cajete 

#endif 
