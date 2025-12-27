#pragma once

#include "parametric_curve.hpp"
#include "bspline_basis.hpp"
#include <array>
#include <vector>

class NURBSCurve : public ParametricCurve
{
public:
    NURBSCurve(std::vector<std::array<double, 3>> control_points,
               std::vector<double> weights,
               std::vector<double> knots,
               int degree,
               bool is_periodic = false);

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
    std::vector<std::array<double, 3>> m_ctrl;
    std::vector<double> m_w;
    std::vector<double> m_knots;
    int m_degree;

    BSplineBasis m_basis;
};