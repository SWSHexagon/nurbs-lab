#include "bspline_surface_builder.hpp"
#include <cmath>

BSplineSurface SurfaceBuilder::Sphere(double R)
{
    // Cubic in both directions
    BSplineBasis ubasis(3, {0, 0, 0, 0, 1, 1, 1, 1});
    BSplineBasis vbasis(3, {0, 0, 0, 0, 1, 1, 1, 1});

    std::vector<std::vector<std::array<double, 3>>> ctrl(
        4,
        std::vector<std::array<double, 3>>(4));

    double c = R * 0.551915024494; // sphere cubic approximation constant

    std::array<std::array<std::array<double, 3>, 4>, 4> P = {{// Row 0
                                                              {{{0, R, 0},
                                                                {c, R, 0},
                                                                {R, c, 0},
                                                                {R, 0, 0}}},
                                                              // Row 1
                                                              {{{0, R, c},
                                                                {c, R, c},
                                                                {R, c, c},
                                                                {R, 0, c}}},
                                                              // Row 2
                                                              {{{0, R, -c},
                                                                {c, R, -c},
                                                                {R, c, -c},
                                                                {R, 0, -c}}},
                                                              // Row 3
                                                              {{{0, R, 0},
                                                                {c, R, 0},
                                                                {R, c, 0},
                                                                {R, 0, 0}}}}};

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ctrl[i][j] = P[i][j];

    return BSplineSurface(ubasis, vbasis, ctrl);
}
