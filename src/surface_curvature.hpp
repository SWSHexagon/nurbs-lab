#pragma once

#include <array>

struct SurfaceCurvature
{
    double K;  // Gaussian
    double H;  // Mean
    double k1; // principal max
    double k2; // principal min
    bool valid;

    std::array<double, 3> dir1; // principal direction for k1 (3D, unit)
    std::array<double, 3> dir2; // principal direction for k2 (3D, unit)

    // optional: UV-space directions
    std::array<double, 2> uv1;
    std::array<double, 2> uv2;
};
