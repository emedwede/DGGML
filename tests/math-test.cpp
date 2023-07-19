#include <iostream>

#include "catch.hpp"

#include "MathUtils.hpp"

TEST_CASE("Running Paramaterized Intersection Test", "[math test]")
{    
    std::size_t N = 3; //number of dimensions

    double x1[3] = {1, 1, 0};
    double x2[3] = {0, 3, 0};
    double x3[3] = {3, 3, 0};
    double u1[3]  = {0, 1, 0}; //unit vector for x1
    double u2[3];

    DGGML::set_unit_vector(x2, x3, u2);

    double sol[2]; //results for intersection parameters in 2D

    double d_p = DGGML::unit_dot_product(u1, u2);

    //ensure the dot product is valid
    REQUIRE(d_p == 0.0);

    DGGML::paramaterized_intersection(x1, x2, x3, u1, sol);
    
    REQUIRE(sol[0] == 2);
    REQUIRE(sol[1] == 1.0/3.0);
    
    //order of x2 and x3 doesn't matter
    DGGML::paramaterized_intersection(x1, x3, x2, u1, sol);
    REQUIRE(sol[0] == 2);
    REQUIRE(sol[1] == 2.0/3.0); //relative intersection point changes
    
    //move the starting point onto the opposite side of x2-x3 and
    //ensure the intersection parameter s[0] is negative
    x1[1] = 5;
    DGGML::paramaterized_intersection(x1, x2, x3, u1, sol);
    REQUIRE(sol[0] == -2);
    REQUIRE(sol[1] == 1.0/3.0);
}
