#pragma once

#include "parametric_curve.hpp"
#include "math_utils.hpp"

using namespace MathUtils;

class LineCurve : public ParametricCurve
{
public:
    LineCurve(const std::array<double, 3> &p0,
              const std::array<double, 3> &p1)
        : p0_(p0), p1_(p1)
    {
        initialize_domain();
    }

    std::array<double, 3> evaluate(double t) const override
    {
        return add(p0_, scale(dir_, t));
    }

    std::array<double, 3> derivative(double t) const override
    {
        return dir_;
    }

    std::array<double, 3> second_derivative(double t) const override
    {
        return {0, 0, 0};
    }

    std::array<double, 3> tangent(double t, bool unitize = true) const override
    {
        return unitize ? normalize(dir_) : dir_;
    }

    std::array<double, 3> normal(double t) const override
    {
        return {0, 0, 0}; // undefined for a line, but safe
    }

    std::pair<double, double> domain() const override
    {
        return {t_min, t_max};
    }

    void DumpInfo() const override
    {
        // print endpoints, direction, domain
    }

protected:
    void initialize_domain() override
    {
        t_min = 0.0;
        t_max = 1.0;
        dir_ = sub(p1_, p0_);
    }

private:
    std::array<double, 3> p0_;
    std::array<double, 3> p1_;
    std::array<double, 3> dir_;
};
