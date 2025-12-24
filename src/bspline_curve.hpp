#pragma once

#include "parametric_curve.hpp"
#include "bspline_basis.hpp"
#include "closest_point_result.hpp"
#include <array>
#include <vector>

class BSplineCurve : public ParametricCurve
{
public:
    BSplineCurve(
        BSplineBasis basis,
        std::vector<std::array<double, 3>> control_points);

    // Basic geometry
    virtual std::array<double, 3> evaluate(double t) const;
    virtual std::array<double, 3> derivative(double t) const;
    virtual std::array<double, 3> second_derivative(double t) const;

    // Convenience
    virtual std::array<double, 3> tangent(double t, bool unitize = true) const;
    virtual std::array<double, 3> normal(double t) const; // principal normal (Frenet)

    // Domain of valid parameters
    virtual std::pair<double, double> domain() const { return {0.0, 1.0}; } // Default is [0, 1]

    double curvature(double t) const;

    virtual ClosestPointResult closest_point_LM(
        const std::array<double, 3> &point,
        double t0,
        int maxIters = 20,
        double tol = 1e-8) const;

    void DumpInfo() const;

    const BSplineBasis &basis() const { return basis_; }
    const std::vector<std::array<double, 3>> &control_points() const { return ctrl_; }

private:
    BSplineBasis basis_;
    std::vector<std::array<double, 3>> ctrl_;

    static constexpr double LAMBDA_MIN = 1e-12;
    static constexpr double LAMBDA_LARGE = 1e6;
    static constexpr double LAMBDA_MAX = 1e8;
    static constexpr double LAMBDA_SEED = 1e-3;
    static constexpr double NEAR_ZERO = 1e-12;

    static void project_to_domain(double &t);

    static double dot(const std::array<double, 3> &a, const std::array<double, 3> &b);
    static std::array<double, 3> sub(const std::array<double, 3> &a, const std::array<double, 3> &b);
    static std::array<double, 3> add(const std::array<double, 3> &a, const std::array<double, 3> &b);
    static std::array<double, 3> scale(const std::array<double, 3> &a, double s);
    static double norm(const std::array<double, 3> &a);
    static std::array<double, 3> normalize(const std::array<double, 3> &a);
};
