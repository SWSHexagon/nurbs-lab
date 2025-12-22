#pragma once

#include <array>

struct SurfaceDifferential
{
    std::array<double, 3> Su;
    std::array<double, 3> Sv;
    std::array<double, 3> Suu;
    std::array<double, 3> Svv;
    std::array<double, 3> Suv;
    std::array<double, 3> N;
    double E;
    double F;
    double G;
    double e;
    double f;
    double g;
    double metricDet;
    bool metricDegenerate;
};
