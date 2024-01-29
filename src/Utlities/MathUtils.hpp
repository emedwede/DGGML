#ifndef DGGML_MATH_UTILS_HPP
#define DGGML_MATH_UTILS_HPP

#include <iostream>
#include <math.h>

namespace DGGML
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

//Order matters here
template <typename PositionType, std::size_t N>
void calculate_difference(PositionType (&p1)[N], PositionType (&p2)[N], PositionType (&p3)[N])
{
    for(auto i = 0; i < N; i++) 
        p3[i] = p1[i] - p2[i];
}

template <typename PositionType, std::size_t N>
double calculate_alpha(PositionType (&pS)[N], PositionType (&pE)[N], PositionType (&pIN)[N])
{
    double result = 0.0;
    for(int i = 0; i < N; i++)
    {
        result += pIN[i]*pE[i] - pS[i] - pS[i]*pE[i] - pS[i];
    }
    result /= calculate_distance(pE, pS);

    return result;
}

template <typename PositionType, std::size_t N>
double calculate_beta(PositionType (&pS)[N], PositionType (&pE)[N], PositionType (&pIN)[N], double alpha)
{
    double result = 0.0;

    double diff = 0.0;

    if(0.0 > alpha && alpha < 1.0)
    {
        for(auto i = 0; i < N; i++)
        {
            double intersection = pS[i]+alpha*(pE[i]-pS[i]);
            diff += pow(pIN[i] - intersection, 2.0);
        }
    } 
    else if(alpha >= 1.0) 
    {
        for(auto i = 0; i < N; i++)
        {
            diff += pow(pIN[i] - pE[i], 2.0);
        }
    }
    else //alpha <= 0.0 
    {
        for(auto i = 0; i < N; i++)
        {
            diff += pow(pIN[i] - pS[i], 2.0);
        }
    }

    return diff;
}


//Order matters here
template<typename PositionType, std::size_t N>
void set_unit_vector(PositionType (&p1)[N], PositionType (&p2)[N], PositionType (&u)[N])
{
    auto len = calculate_distance(p1, p2);
    calculate_difference(p1, p2, u);
    for(auto i = 0; i < N; i++) 
        u[i] /= len;
}

template<typename PositionType, std::size_t N>
void paramaterized_intersection(
        PositionType (&p1)[N], 
        PositionType (&p2)[N],
        PositionType (&p3)[N],
        PositionType (&u)[N], 
        PositionType (&s)[2])
{
    //Calculation copied from mathematica, TODO: make more efficient 
    s[0] = (-p2[1]*p3[0]+p1[1]*(-p2[0]+p3[0])+p1[0]*(p2[1]-p3[1])+p2[0]*p3[1])
            /
            (-p2[1]*u[0]+p3[1]*u[0]+(p2[0]-p3[0])*u[1]);
    
    s[1] = (p1[1]*u[0]-p2[1]*u[0]+(-p1[0]+p2[0])*u[1])
            /
            (-p2[1]*u[0]+p3[1]*u[0]+(p2[0]-p3[0])*u[1]);

}


double cross_product(double a1, double a2, double b1, double b2)
{
    return ( ( a1 * b2 ) - ( a2 * b1 ) );
}

template <typename UnitVecType, std::size_t N>
double unit_dot_product(UnitVecType (&u1)[N], UnitVecType (&u2)[N])
{
    double result = 0.0;

    for(auto i = 0; i < N; i++) {
        result += u1[i]*u2[i];
    }

    return result;
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

//calculates the distance of point Q to a line segment between two points P1 and P2
double distanceToLineSegment(double x1, double y1, double x2, double y2, double xq, double yq) {

    // Ensure x1 <= x2, swap if needed
    if (x1 > x2) {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    // Calculate parametric value (t)
    double t = ((xq - x1) * (x2 - x1) + (yq - y1) * (y2 - y1)) / ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Check if the projection is within the line segment
    double xc, yc;
    if (t < 0) {
        xc = x1;
        yc = y1;
    } else if (t > 1) {
        xc = x2;
        yc = y2;
    } else {
        // Calculate closest point on the line
        xc = x1 + t * (x2 - x1);
        yc = y1 + t * (y2 - y1);
    }

    // Calculate distance between the point and the closest point on the line
    double distance = sqrt((xq - xc) * (xq - xc) + (yq - yc) * (yq - yc));
    return distance;
}

} //end namespace DGGML

#endif 
