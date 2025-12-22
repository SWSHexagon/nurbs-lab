#pragma once

#include "bspline_surface.hpp"

namespace SurfaceBuilder
{
    BSplineSurface Bowl();
    BSplineSurface Sphere(double radius);
    BSplineSurface Cylinder(double radius, double height);
    BSplineSurface Plane(
        double width, double height,
        const std::array<double, 3> &origin,
        const std::array<double, 3> &udir,
        const std::array<double, 3> &vdir);
}