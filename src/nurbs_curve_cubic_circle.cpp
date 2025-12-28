#include "nurbs_curve_builder.hpp"
#include "math_utils.hpp"
#include <numbers>
#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>

using namespace MathUtils;

NURBSCurve NURBSCurveBuilder::MakeNURBSCubicCircle(
    const std::array<double, 3> &center,
    const std::array<double, 3> &normal,
    double radius)
{
    using Vec3 = std::array<double, 3>;

#if 0
    // Build orthonormal frame (u, v, n)
    auto n = normalize(normal);
    auto tmp = (std::abs(n[0]) < 0.9) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    auto u = normalize(cross(n, tmp));
    auto v = cross(n, u);
#endif

    // Build control points and weights
    std::vector<Vec3> ctrl =
        {
            {1.035276, 0.000000, 0.000000},
            {0.896575, 0.517638, 0.000000},
            {0.517638, 0.896575, 0.000000},
            {0.000000, 1.035276, 0.000000},
            {-0.517638, 0.896575, 0.000000},
            {-0.896575, 0.517638, 0.000000},
            {-1.035276, 0.000000, 0.000000},
            {-0.896575, -0.517638, 0.000000},
            {-0.517638, -0.896575, 0.000000},
            {0.000000, -1.035276, 0.000000},
            {0.517638, -0.896575, 0.000000},
            {0.896575, -0.517638, 0.000000},
            {1.035276, 0.000000, 0.000000},
            {0.896575, 0.517638, 0.000000},
            {0.517638, 0.896575, 0.000000}};

    std::vector<double> w =
        {
            1.000000,
            0.965926,
            0.965926,
            1.000000,
            0.965926,
            0.965926,
            1.000000,
            0.965926,
            0.965926,
            1.000000,
            0.965926,
            0.965926,
            1.000000,
            0.965926,
            0.965926};

    // Correct periodic cubic knot vector
    std::vector<double> knots =
        {
            -0.25,
            -0.1666667,
            -0.0833333,
            0.0000000,
            0.0833333,
            0.1666667,
            0.2500000,
            0.3333333,
            0.4166667,
            0.5000000,
            0.5833333,
            0.6666667,
            0.7500000,
            0.8333333,
            0.9166667,
            1.0000000,
            1.0833333,
            1.1666667,
            1.2500000};

    const int degree = 3;
    bool isPeriodic = true;

    return NURBSCurve(ctrl, w, knots, degree, isPeriodic);
}
