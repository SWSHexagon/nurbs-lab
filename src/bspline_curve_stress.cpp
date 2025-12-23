#include "bspline_curve_builder.hpp"

BSplineCurve CurveBuilder::Stress()
{
    const int degree = 3;
    const int numCtrl = 10;

    // 1. Control points: alternating high curvature
    // Control points: moderate oscillation + gentle 3D twist
    std::vector<std::array<double, 3>> ctrl =
        {
            {0.0, 0.0, 0.0},
            {1.0, 1.5, 0.5},
            {2.0, -1.5, 1.0},
            {3.0, 1.2, 1.5},
            {4.0, -1.2, 2.0},
            {5.0, 1.0, 2.5},
            {6.0, -1.0, 3.0},
            {7.0, 0.8, 3.5},
            {8.0, -0.8, 4.0},
            {9.0, 0.0, 4.5}};

    // 2. Knot vector: open-uniform, degree 3
    //    Total knots = numCtrl + degree + 1 = 10 + 3 + 1 = 14
    std::vector<double> knots(14);

    // First 4 knots = 0
    for (int i = 0; i < 4; ++i)
        knots[i] = 0.0;

    // Internal knots: 6 evenly spaced values
    // internal count = numCtrl - degree - 1 = 10 - 3 - 1 = 6
    for (int i = 1; i <= 6; ++i)
        knots[3 + i] = double(i) / 7.0; // 1/7, 2/7, ..., 6/7

    // Last 4 knots = 1
    for (int i = 0; i < 4; ++i)
        knots[10 + i] = 1.0;

    // 3. Construct the curve
    BSplineBasis basis(degree, knots);
    return BSplineCurve(basis, ctrl);
}
