#pragma once

#include <array>
#include <utility>
#include "closest_point_result.hpp"

class ParametricCurve
{
public:
    virtual ~ParametricCurve() {}

    // Core evaluation
    virtual std::array<double, 3> evaluate(double t) const = 0;
    virtual std::array<double, 3> derivative(double t) const = 0;
    virtual std::array<double, 3> second_derivative(double t) const = 0;

    // Convenience
    virtual std::array<double, 3> tangent(double t, bool unitize = true) const = 0; // unit tangent (Frenet)
    virtual std::array<double, 3> normal(double t) const = 0;                       // principal normal (Frenet)

    // Domain of valid parameters
    virtual std::pair<double, double> domain() const = 0;

    virtual ClosestPointResult closest_point_LM(
        const std::array<double, 3> &point,
        double t0,
        int maxIters = 20,
        double tol = 1e-8) const = 0;
};