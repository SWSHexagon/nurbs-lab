#include "bspline_surface_builder.hpp"

BSplineSurface SurfaceBuilder::Bowl()
{
    BSplineBasis ub(2, {0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1});
    BSplineBasis vb(2, {0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1});

    // 4×4 tilted curved control net
    std::vector<std::vector<std::array<double, 3>>> ctrl(5, std::vector<std::array<double, 3>>(5));

    const double ax = 0;   // tilt in x
    const double ay = 0;   // tilt in y
    const double k = 0.25; // curvature strength

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            double x = double(i) - 2.0;
            double y = double(j) - 2.0;
            double z = k * (x * x + y * y) + ax * x + ay * y;
            ctrl[i][j] = {x, y, z};
        }
    }

    return BSplineSurface(ub, vb, ctrl);
}
