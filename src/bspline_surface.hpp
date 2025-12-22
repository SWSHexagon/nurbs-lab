#pragma once
#include "bspline_basis.hpp"
#include <vector>
#include <array>
#include <tuple>

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
    std::array<double, 3> normal(double u, double v, bool unitize = true) const;

    struct ClosestPointResult
    {
        double u = 0.0;
        double v = 0.0;

        std::array<double, 3> point3D; // final closest point on surface

        double distance = 0.0; // final |S(u,v) - P|
        double gradNorm = 0.0; // ||∇(0.5||F||²)||

        int iterations = 0; // final number of iterations

        enum class Status
        {
            Success,      // Interior convergence
            Boundary,     // Converged on boundary
            Stagnation,   // Step size too small with large gradient norm
            Divergence,   // lambda too large or Hessian unusable
            MaxIterations // Loop exhaused without convergence
        } status = Status::MaxIterations;

        bool onBoundary = false; // Whether the final point is on the boundary
    };

    // Returns the closest point on the surface to the given point in space
    // using Levenberg-Marquardt gradient descent.
    ClosestPointResult closest_point_LM(
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
        double metricDet;
        bool metricDegenerate;
    };

    SurfaceDifferential differential(double u, double v) const;

    struct Curvature
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

    Curvature curvature(double u, double v) const;

    void DumpInfo() const;

private:
    BSplineBasis u_basis_;
    BSplineBasis v_basis_;
    std::vector<std::vector<std::array<double, 3>>> ctrl_;

    static const double LAMBDA_MIN;
    static const double LAMBDA_LARGE;
    static const double LAMBDA_MAX;
    static const double LAMBDA_SEED;
    static const double NEAR_ZERO;

    static void project_to_domain(double &u, double &v);

    static std::array<double, 2> principal_direction_uv(
        const SurfaceDifferential &d,
        double k);

    static std::array<double, 3> combine_uv_to_3d(
        const SurfaceDifferential &d,
        const std::array<double, 2> &uv);
};