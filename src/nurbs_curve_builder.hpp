#pragma once

#include "nurbs_curve.hpp"

namespace NURBSCurveBuilder
{
    NURBSCurve MakeNURBSCircle(
        const std::array<double, 3> &center,
        const std::array<double, 3> &normal,
        double radius);

    NURBSCurve MakeNURBSPeriodicCircle(double radius);

    NURBSCurve MakeNURBSCubicCircle(
        const std::array<double, 3> &center,
        const std::array<double, 3> &normal,
        double radius);
}