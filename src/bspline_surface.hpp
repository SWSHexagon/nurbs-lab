#pragma once
#include "bspline_basis.hpp"
#include <vector>
#include <array>

class BSplineSurface
{
public:
    BSplineSurface(
        BSplineBasis u_basis,
        BSplineBasis v_basis,
        std::vector<std::vector<std::array<double, 3>>> control_points);

    std::array<double, 3> evaluate(double u, double v) const;
    std::array<double, 3> derivative_u(double u, double v) const;
    std::array<double, 3> derivative_v(double u, double v) const;
    std::array<double, 3> second_derivative_uu(double u, double v) const;
    std::array<double, 3> second_derivative_vv(double u, double v) const;
    std::array<double, 3> second_derivative_uv(double u, double v) const;
    std::array<double, 3> normal(double u, double v) const;

    std::pair<double, double> closest_point(
        const std::array<double, 3> &point,
        double u0,
        double v0,
        int maxIters = 20,
        double tol = 1e-6) const;

    std::pair<double, double> closest_point_LM(
        const std::array<double, 3> &point,
        double u0,
        double v0,
        int maxIters = 20,
        double tol = 1e-6) const;

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
    };

    SurfaceDifferential differential(double u, double v) const;

    struct Curvature
    {
        double K;  // Gaussian
        double H;  // Mean
        double k1; // principal max
        double k2; // principal min
        bool valid;
    };

    Curvature curvature(double u, double v) const;

    void DumpInfo() const;

private:
    BSplineBasis u_basis_;
    BSplineBasis v_basis_;
    std::vector<std::vector<std::array<double, 3>>> ctrl_;
};