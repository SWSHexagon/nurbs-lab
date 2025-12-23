#pragma once

#include "bspline_basis.hpp"
#include "surface_differential.hpp"
#include "surface_curvature.hpp"
#include "closest_point_result.hpp"
#include "surface_frame.hpp"
#include "bspline_curve.hpp"
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

    // Returns the closest point on the surface to the given point in space
    // using Levenberg-Marquardt gradient descent.
    ClosestPointResult closest_point_LM(
        const std::array<double, 3> &point,
        double u0 = 0.5,
        double v0 = 0.5,
        int maxIters = 20,
        double tol = 1e-8) const;

    ClosestPointResult closest_point_with_boundary_fallback(
        const std::array<double, 3> &point,
        double u0 = 0.5,
        double v0 = 0.5,
        int maxIters = 20,
        double tol = 1e-8,
        int maxCurveIter = 20,
        double curveTol = 1e-8) const;

    ClosestPointResult closest_point_global(
        const std::array<double, 3> &point,
        int gridResolution = 3, // e.g. 3 or 4
        double u0 = 0.5,
        double v0 = 0.5,
        int maxIter = 20,
        double tol = 1e-8,
        int maxCurveIter = 20,
        double curveTol = 1e-8) const;

    SurfaceDifferential differential(double u, double v) const;

    SurfaceCurvature curvature(double u, double v) const;

    SurfaceFrame frame(double u, double v) const;

    const BSplineCurve &boundary_u_min() const { return boundary_u_min_; }
    const BSplineCurve &boundary_u_max() const { return boundary_u_max_; }
    const BSplineCurve &boundary_v_min() const { return boundary_v_min_; }
    const BSplineCurve &boundary_v_max() const { return boundary_v_max_; }

    void DumpInfo() const;

private:
    BSplineBasis u_basis_;
    BSplineBasis v_basis_;
    std::vector<std::vector<std::array<double, 3>>> ctrl_;

    BSplineCurve boundary_u_min_{BSplineBasis(0, {0}), {}};
    BSplineCurve boundary_u_max_{BSplineBasis(0, {0}), {}};
    BSplineCurve boundary_v_min_{BSplineBasis(0, {0}), {}};
    BSplineCurve boundary_v_max_{BSplineBasis(0, {0}), {}};

    void build_boundary_curves();

    static constexpr double LAMBDA_MIN = 1e-12;
    static constexpr double LAMBDA_LARGE = 1e6;
    static constexpr double LAMBDA_MAX = 1e8;
    static constexpr double LAMBDA_SEED = 1e-3;
    static constexpr double NEAR_ZERO = 1e-12;

    static void project_to_domain(double &u, double &v);

    static std::array<double, 2> principal_direction_uv(
        const SurfaceDifferential &d,
        double k);

    static std::array<double, 3> combine_uv_to_3d(
        const SurfaceDifferential &d,
        const std::array<double, 2> &uv);
};