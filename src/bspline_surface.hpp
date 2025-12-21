#pragma once
#include "bspline_basis.hpp"
#include <vector>
#include <array>

class BSplineSurface {
public:
    BSplineSurface(
        BSplineBasis u_basis,
        BSplineBasis v_basis,
        std::vector<std::vector<std::array<double,3>>> control_points
    );

    std::array<double,3> evaluate(double u, double v) const;

private:
    BSplineBasis u_basis_;
    BSplineBasis v_basis_;
    std::vector<std::vector<std::array<double,3>>> ctrl_;
};