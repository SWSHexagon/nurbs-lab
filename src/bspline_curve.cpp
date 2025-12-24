#include "bspline_curve.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

BSplineCurve::BSplineCurve(
    BSplineBasis basis,
    std::vector<std::array<double, 3>> control_points)
    : basis_(std::move(basis)), ctrl_(std::move(control_points))
{
}

void BSplineCurve::initialize_domain()
{
    // Initialize curve domain
    auto [tMin, tMax] = basis_.getExtents();

    t_min = tMin;
    t_max = tMax;
}

std::array<double, 3> BSplineCurve::evaluate(double t) const
{
    std::array<double, 3> C{0, 0, 0};

    int n = static_cast<int>(ctrl_.size());
    int p = basis_.degree();

    for (int i = 0; i < n; ++i)
    {
        double Ni = basis_.evaluate(i, p, t);
        C[0] += Ni * ctrl_[i][0];
        C[1] += Ni * ctrl_[i][1];
        C[2] += Ni * ctrl_[i][2];
    }

    return C;
}

std::array<double, 3> BSplineCurve::derivative(double t) const
{
    std::array<double, 3> Cp{0, 0, 0};

    int n = basis_.numBasis();
    for (int i = 0; i < n; ++i)
    {
        double Ni_t = basis_.derivative(i, t);

        Cp[0] += Ni_t * ctrl_[i][0];
        Cp[1] += Ni_t * ctrl_[i][1];
        Cp[2] += Ni_t * ctrl_[i][2];
    }

    return Cp;
}

std::array<double, 3> BSplineCurve::second_derivative(double t) const
{
    std::array<double, 3> Cpp{0, 0, 0};

    int n = basis_.numBasis();
    for (int i = 0; i < n; ++i)
    {
        double Ni_tt = basis_.second_derivative(i, t);

        Cpp[0] += Ni_tt * ctrl_[i][0];
        Cpp[1] += Ni_tt * ctrl_[i][1];
        Cpp[2] += Ni_tt * ctrl_[i][2];
    }

    return Cpp;
}

std::array<double, 3> BSplineCurve::tangent(double t, bool unitize) const
{
    auto Cp = derivative(t);
    if (unitize)
        return normalize(Cp);
    return Cp;
}

double BSplineCurve::curvature(double t) const
{
    auto Cp = derivative(t);
    auto Cpp = second_derivative(t);

    double speed2 = dot(Cp, Cp);
    double speed = std::sqrt(speed2);

    if (speed < NEAR_ZERO)
        return 0.0;

    // kappa = |C' x C''| / |C'|^3
    std::array<double, 3> cross{
        Cp[1] * Cpp[2] - Cp[2] * Cpp[1],
        Cp[2] * Cpp[0] - Cp[0] * Cpp[2],
        Cp[0] * Cpp[1] - Cp[1] * Cpp[0]};

    double num = norm(cross);
    double denom = speed2 * speed; // |C'|^3

    if (denom < NEAR_ZERO)
        return 0.0;

    return num / denom;
}

std::array<double, 3> BSplineCurve::normal(double t) const
{
    auto Cp = derivative(t);
    auto Cpp = second_derivative(t);

    double speed2 = dot(Cp, Cp);
    if (speed2 < NEAR_ZERO)
        return {0, 0, 0};

    // N = (C'' - (C''·T)T) normalized
    auto T = normalize(Cp);
    double proj = dot(Cpp, T);

    auto N = sub(Cpp, scale(T, proj));
    return normalize(N);
}

void BSplineCurve::DumpInfo() const
{
    std::cout << "BSplineCurve:\n";
    std::cout << "  Degree: " << basis_.degree() << "\n";
    std::cout << "  Num ctrl pts: " << ctrl_.size() << "\n";
    std::cout << "  Knots: ";
    for (double k : basis_.knots())
        std::cout << k << " ";
    std::cout << "\n";

    std::cout << "  Control points:\n";
    for (size_t i = 0; i < ctrl_.size(); ++i)
    {
        std::cout << "    P[" << i << "] = ("
                  << ctrl_[i][0] << ", "
                  << ctrl_[i][1] << ", "
                  << ctrl_[i][2] << ")\n";
    }
}