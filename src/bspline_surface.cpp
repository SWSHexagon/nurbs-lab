#include "bspline_surface.hpp"

BSplineSurface::BSplineSurface(
    BSplineBasis u_basis,
    BSplineBasis v_basis,
    std::vector<std::vector<std::array<double, 3>>> control_points)
    : u_basis_(std::move(u_basis)),
      v_basis_(std::move(v_basis)),
      ctrl_(std::move(control_points))
{
}

std::array<double, 3> BSplineSurface::evaluate(double u, double v) const
{
    std::array<double, 3> p = {0.0, 0.0, 0.0};

    int nu = static_cast<int>(ctrl_.size());
    int nv = static_cast<int>(ctrl_[0].size());

    for (int i = 0; i < nu; ++i)
    {
        double Nu = u_basis_.evaluate(i, u_basis_.degree(), u);

        for (int j = 0; j < nv; ++j)
        {
            double Nv = v_basis_.evaluate(j, v_basis_.degree(), v);

            double w = Nu * Nv;

            p[0] += w * ctrl_[i][j][0];
            p[1] += w * ctrl_[i][j][1];
            p[2] += w * ctrl_[i][j][2];
        }
    }

    return p;
}
