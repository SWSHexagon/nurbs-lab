#include "bspline_surface_builder.hpp"

BSplineSurface SurfaceBuilder::Plane(
    double width, double height,
    const std::array<double, 3> &origin,
    const std::array<double, 3> &udir,
    const std::array<double, 3> &vdir)
{
    BSplineBasis ubasis(1, {0, 0, 1, 1});
    BSplineBasis vbasis(1, {0, 0, 1, 1});

    std::vector<std::vector<std::array<double, 3>>> ctrl(2,
                                                         std::vector<std::array<double, 3>>(2));

    ctrl[0][0] = origin;
    ctrl[1][0] = {origin[0] + width * udir[0],
                  origin[1] + width * udir[1],
                  origin[2] + width * udir[2]};

    ctrl[0][1] = {origin[0] + height * vdir[0],
                  origin[1] + height * vdir[1],
                  origin[2] + height * vdir[2]};

    ctrl[1][1] = {origin[0] + width * udir[0] + height * vdir[0],
                  origin[1] + width * udir[1] + height * vdir[1],
                  origin[2] + width * udir[2] + height * vdir[2]};

    return BSplineSurface(ubasis, vbasis, ctrl);
}
