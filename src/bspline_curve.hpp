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
    virtual std::pair<double, double> domain() const { return {t_min, t_max}; } // Default is [0, 1]

    double curvature(double t) const;

    void DumpInfo() const;

    const BSplineBasis &basis() const { return basis_; }
    const std::vector<std::array<double, 3>> &control_points() const { return ctrl_; }

protected:
    void initialize_domain();

private:
    BSplineBasis basis_;
    std::vector<std::array<double, 3>> ctrl_;
};
