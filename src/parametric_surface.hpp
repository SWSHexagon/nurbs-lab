#pragma once

#include "closest_point_result.hpp"
#include "parametric_curve.hpp"

#include <array>
#include <utility>

class ParametricSurface
{
public:
    virtual ~ParametricSurface() {}

    // Core evaluation
    virtual std::array<double, 3> evaluate(double u, double v) const = 0;

    // First partials
    virtual std::pair<std::array<double, 3>, std::array<double, 3>> derivatives(
        double u, double v) const = 0;

    // Normal vector
    virtual std::array<double, 3> normal(double u, double v, bool unitize = true) const = 0;

    // Domain rectangle
    virtual std::array<double, 4> domain() const = 0; // {u0, u1, v0, v1}

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

    virtual const ParametricCurve &boundary_u_min() const = 0;
    virtual const ParametricCurve &boundary_u_max() const = 0;
    virtual const ParametricCurve &boundary_v_min() const = 0;
    virtual const ParametricCurve &boundary_v_max() const = 0;

protected:
    static constexpr double LAMBDA_MIN = 1e-12;
    static constexpr double LAMBDA_LARGE = 1e6;
    static constexpr double LAMBDA_MAX = 1e8;
    static constexpr double LAMBDA_SEED = 1e-3;
    static constexpr double NEAR_ZERO = 1e-12;

    double u_min = 0;
    double u_max = 1;
    double v_min = 0;
    double v_max = 1;

    void project_to_domain(double &u, double &v) const;
    virtual void initialize_domain() = 0;
};