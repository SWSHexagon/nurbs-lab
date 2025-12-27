#include "nurbs_curve_builder.hpp"
#include "math_utils.hpp"
#include <numbers>
#include <cmath>
#include <array>
#include <iostream>

using namespace MathUtils;

NURBSCurve NURBSCurveBuilder::MakeNURBSPeriodicCircle(double radius)
{
    const double c = std::sqrt(0.5);

    std::vector<std::array<double, 3>> P = {
        {radius, 0, 0},
        {radius * c, radius * c, 0},
        {0, radius, 0},
        {-radius * c, radius * c, 0},
        {-radius, 0, 0},
        {-radius * c, -radius * c, 0},
        {0, -radius, 0},
        {radius * c, -radius * c, 0},
        {radius, 0, 0} // wrap
    };

    std::vector<double> W = {
        1, c, 1, c, 1, c, 1, c, 1};

    std::vector<double> U(12);
    for (int i = 0; i < 12; ++i)
        U[i] = static_cast<double>(i);

    return NURBSCurve(P, W, U, 2, true);
}
