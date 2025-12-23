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

double BSplineCurve::dot(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> BSplineCurve::sub(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::array<double, 3> BSplineCurve::add(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::array<double, 3> BSplineCurve::scale(const std::array<double, 3> &a, double s)
{
    return {a[0] * s, a[1] * s, a[2] * s};
}

double BSplineCurve::norm(const std::array<double, 3> &a)
{
    return std::sqrt(dot(a, a));
}

std::array<double, 3> BSplineCurve::normalize(const std::array<double, 3> &a)
{
    double n = norm(a);
    if (n < NEAR_ZERO)
        return {0, 0, 0};
    return {a[0] / n, a[1] / n, a[2] / n};
}

void BSplineCurve::project_to_domain(double &t)
{
    // Assuming basis domain is [0,1] mapped from knots
    if (t < 0.0)
        t = 0.0;
    if (t > 1.0)
        t = 1.0;
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

BSplineCurve::ClosestPointResult BSplineCurve::closest_point_LM(
    const std::array<double, 3> &P,
    double t0,
    int maxIters,
    double tol) const
{
    ClosestPointResult result;
    result.t = t0;

    double lambda = LAMBDA_SEED;

    for (int iter = 0; iter < maxIters; ++iter)
    {
        project_to_domain(result.t);

        auto C = evaluate(result.t);
        auto Cp = derivative(result.t);
        auto Cpp = second_derivative(result.t);

        auto r = sub(C, P);         // residual
        double f = 0.5 * dot(r, r); // objective

        double g = dot(r, Cp);                // d/dt(0.5|r|^2) = r·C'
        double H = dot(Cp, Cp) + dot(r, Cpp); // scalar Hessian

        result.gradNorm = std::abs(g);

        // Convergence check
        if (result.gradNorm < tol)
        {
            result.point3D = C;
            result.distance = std::sqrt(2.0 * f);
            result.iterations = iter;
            result.status = ClosestPointResult::Status::Success;
            result.onBoundary = (result.t <= 0.0 + 1e-10 || result.t >= 1.0 - 1e-10);
            return result;
        }

        // Damped scalar step: (H + lambda) * delta = -g
        double denom = H + lambda;
        if (std::abs(denom) < NEAR_ZERO)
        {
            // Hessian unusable: treat as divergence
            result.status = ClosestPointResult::Status::Divergence;
            break;
        }

        double delta = -g / denom;

        if (std::abs(delta) < tol * 1e-2)
        {
            // Stagnation
            result.status = ClosestPointResult::Status::Stagnation;
            break;
        }

        double t_new = result.t + delta;
        project_to_domain(t_new);

        auto C_new = evaluate(t_new);
        auto r_new = sub(C_new, P);
        double fNew = 0.5 * dot(r_new, r_new);

        // Check improvement
        if (fNew < f)
        {
            // Accept step, decrease lambda
            result.t = t_new;
            lambda /= 10;
            lambda = std::max(lambda, LAMBDA_MIN);
        }
        else
        {
            // Reject step, increase lambda
            lambda *= 10;
            lambda = std::min(lambda, LAMBDA_MAX);
        }

        result.iterations = iter;
    }

    // Final evaluation
    project_to_domain(result.t);
    auto C = evaluate(result.t);
    auto r = sub(C, P);

    result.point3D = C;
    result.distance = norm(r);
    result.onBoundary = (result.t <= 0.0 + 1e-10 || result.t >= 1.0 - 1e-10);

    if (result.status == ClosestPointResult::Status::MaxIterations)
        result.status = ClosestPointResult::Status::MaxIterations;

    return result;
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