#pragma once

#include <array>

struct SurfaceFrame
{
    // Position on the surface
    std::array<double, 3> P;

    // Tangent vectors
    std::array<double, 3> Su;
    std::array<double, 3> Sv;

    // Normal (unit)
    std::array<double, 3> N;

    // Principal curvatures
    double k1;
    double k2;

    // Principal directions (3D, unit)
    std::array<double, 3> dir1;
    std::array<double, 3> dir2;

    // Optional: UV-space principal directions
    std::array<double, 2> uv1;
    std::array<double, 2> uv2;

    // First fundamental form
    double E, F, G;

    // Second fundamental form
    double e, f, g;

    // Metric determinant and degeneracy flag
    double metricDet;
    bool metricDegenerate;

    // Curvature validity
    bool valid;
};
