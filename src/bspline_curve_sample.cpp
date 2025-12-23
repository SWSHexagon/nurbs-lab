#include "bspline_curve_builder.hpp"
#include <cmath>
#include <vector>
#include <array>

BSplineCurve CurveBuilder::Sample()
{
    const int degree = 3;
    const int numCtrl = 10;

    // 1. Control points (example: random or structured)
    std::vector<std::array<double, 3>> ctrl(numCtrl);

    for (int i = 0; i < numCtrl; ++i)
    {
        double x = double(i);
        double y = std::sin(i * 0.4);
        double z = 0.1 * i;
        ctrl[i] = {x, y, z};
    }

    // 2. Knot vector (open-uniform)
    std::vector<double> knots;
    knots.resize(numCtrl + degree + 1); // 10 + 3 + 1 = 14

    // First 4 knots = 0
    for (int i = 0; i < degree + 1; ++i)
        knots[i] = 0.0;

    // Internal knots
    int internal = numCtrl - degree - 1; // 10 - 3 - 1 = 6
    for (int i = 1; i <= internal; ++i)
        knots[degree + i] = double(i) / (internal + 1);

    // Last 4 knots = 1
    for (int i = 0; i < degree + 1; ++i)
        knots[numCtrl + i] = 1.0;

    // 3. Build curve
    BSplineBasis basis(degree, knots);
    return BSplineCurve(basis, ctrl);
}
