#pragma once

#include "closest_point_result.hpp"
#include "math_utils.hpp"
#include <array>
#include <utility>

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
        double tol = 1e-8) const;

    virtual void DumpInfo() const = 0;

protected:
    static constexpr double LAMBDA_MIN = 1e-12;
    static constexpr double LAMBDA_LARGE = 1e6;
    static constexpr double LAMBDA_MAX = 1e8;
    static constexpr double LAMBDA_SEED = 1e-3;
    static constexpr double NEAR_ZERO = 1e-12;

    double t_min = 0;
    double t_max = 1;

    bool m_is_periodic = false;

    void project_to_domain(double &t) const;
    virtual void initialize_domain() = 0;

    void set_periodic(bool is_Periodic) { m_is_periodic = is_Periodic; }
};
