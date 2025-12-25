#include "nurbs_curve_builder.hpp"
#include "math_utils.hpp"
#include <numbers>
#include <cmath>
#include <array>

using namespace MathUtils;

NURBSCurve NURBSCurveBuilder::MakeNURBSCircle(
    const std::array<double, 3> &center,
    const std::array<double, 3> &normal,
    double radius)
{
    const double pi = std::numbers::pi;

    // Build orthonormal frame
    auto n = normalize(normal);
    auto tmp = (std::abs(n[0]) < 0.9)
                   ? std::array<double, 3>{1, 0, 0}
                   : std::array<double, 3>{0, 1, 0};

    auto u = normalize(cross(n, tmp));
    auto v = cross(n, u);

    // 8 control points (no duplication)
    std::vector<std::array<double, 3>> ctrl(8);
    std::vector<double> w(8);

    const double w_diag = std::sqrt(0.5);
    const double dtheta = pi / 4.0;

    for (int i = 0; i < 8; ++i)
    {
        double theta = i * dtheta;
        double c = std::cos(theta);
        double s = std::sin(theta);

        ctrl[i] = add(center,
                      add(scale(u, radius * c),
                          scale(v, radius * s)));

        w[i] = (i % 2 == 0) ? 1.0 : w_diag;
    }

    // Clamped knot vector for degree 2, 8 control points
    std::vector<double> knots = {
        0.0, 0.0, 0.0,
        0.25, 0.25,
        0.5, 0.5,
        0.75, 0.75,
        1.0, 1.0, 1.0};

    return NURBSCurve(ctrl, w, knots, 2);
}
