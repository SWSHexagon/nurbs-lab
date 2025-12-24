#pragma once

#include "parametric_curve.hpp"

#pragma once
#include "parametric_curve.hpp"
#include <array>
#include <cmath>

class CircleCurve : public ParametricCurve
{
public:
    CircleCurve(const std::array<double, 3> &center,
                const std::array<double, 3> &normal,
                double radius);

    std::array<double, 3> evaluate(double t) const override;
    std::array<double, 3> derivative(double t) const override;
    std::array<double, 3> second_derivative(double t) const override;

    std::array<double, 3> tangent(double t, bool unitize = true) const override;
    std::array<double, 3> normal(double t) const override;

    std::pair<double, double> domain() const override
    {
        return {t_min, t_max};
    }

    void DumpInfo() const override;

protected:
    void initialize_domain() override;

private:
    std::array<double, 3> C; // center
    std::array<double, 3> N; // unit normal
    std::array<double, 3> U; // orthonormal basis vector
    std::array<double, 3> V; // orthonormal basis vector
    double R;                // radius
};