#include "nurbs_curve_builder.hpp"
#include "math_utils.hpp"
#include <numbers>
#include <cmath>
#include <array>
#include <iostream>

using namespace MathUtils;

NURBSCurve NURBSCurveBuilder::MakeNURBSCircle(
    const std::array<double, 3> &center,
    const std::array<double, 3> &normal,
    double radius)
{
    const double pi = std::numbers::pi;

    // Build orthonormal frame
    auto n = normalize(normal);
    auto tmp = (std::abs(n[0]) < 0.9)
                   ? std::array<double, 3>{1, 0, 0}
                   : std::array<double, 3>{0, 1, 0};

    auto u = normalize(cross(n, tmp));
    auto v = cross(n, u);

    // 8 control points (no duplication)
    std::vector<std::array<double, 3>> ctrl = {
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
        {-1, 1, 0},
        {-1, 0, 0},
        {-1, -1, 0},
        {0, -1, 0},
        {1, -1, 0},
        {1, 0, 0}};

    std::vector<double> w(9);
    const double w_diag = std::sqrt(2) / 2;

    for (int i = 0; i < 9; ++i)
    {
        // ctrl[i] = add(center, scale(ctrl[i], radius));
        auto px = ctrl[i][0];
        auto py = ctrl[i][1];

        ctrl[i] = add(center,
                      add(scale(u, radius * px),
                          scale(v, radius * py)));

        w[i] = (i % 2 == 0) ? 1.0 : w_diag;
    }

    // Clamped knot vector for degree 2, 8 control points
    std::vector<double> knots = {
        0.0, 0.0, 0.0,
        0.25, 0.25,
        0.5, 0.5,
        0.75, 0.75,
        1.0, 1.0, 1.0};

    std::cout << "Knots: ";
    for (const auto &knot : knots)
        std::cout << knot << " ";
    std::cout << std::endl;

    std::cout << "Weights: ";
    for (const auto &weight : w)
        std::cout << weight << " ";
    std::cout << std::endl;

    std::cout << "Control: ";
    for (const auto &p : ctrl)
        std::cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ") ";
    std::cout << std::endl;

    std::cout << "Degree: 2" << std::endl;

    return NURBSCurve(ctrl, w, knots, 2, true);
}
