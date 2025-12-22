#include "bspline_surface_builder.hpp"
#include <cmath>

BSplineSurface SurfaceBuilder::Cylinder(double R, double height)
{
    BSplineBasis ubasis(3, {0, 0, 0, 0, 1, 1, 1, 1});
    BSplineBasis vbasis(1, {0, 0, 1, 1});

    std::vector<std::vector<std::array<double, 3>>> ctrl(4,
                                                         std::vector<std::array<double, 3>>(2));

    double a = R;
    double b = R * std::sqrt(2) / 2;

    // Bottom ring
    ctrl[0][0] = {a, 0, 0};
    ctrl[1][0] = {b, b, 0};
    ctrl[2][0] = {0, a, 0};
    ctrl[3][0] = {-b, b, 0};

    // Top ring
    ctrl[0][1] = {a, 0, height};
    ctrl[1][1] = {b, b, height};
    ctrl[2][1] = {0, a, height};
    ctrl[3][1] = {-b, b, height};

    return BSplineSurface(ubasis, vbasis, ctrl);
}
