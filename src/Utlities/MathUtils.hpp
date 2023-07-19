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

} //end namespace DGGML

#endif 
