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

ClosestPointResult BSplineCurve::closest_point_LM(
    const std::array<double, 3> &P,
    double t0,
    int maxIters,
    double tol) const
{
    ClosestPointResult result{};

    // Current parameter
    double t = t0;
    project_to_domain(t);

    // Evaluate at initial guess
    auto C0 = evaluate(t);
    auto r0 = sub(C0, P);
    double f_best = 0.5 * dot(r0, r0);

    double t_best = t;
    std::array<double, 3> C_best = C0;

    // Seed lambda
    double lambda = LAMBDA_SEED;

    int iter = 0;

    for (iter = 0; iter < maxIters; ++iter)
    {
        project_to_domain(t);

        // Curve and derivatives at current t
        auto C = evaluate(t);
        auto Cp = derivative(t);
        auto Cpp = second_derivative(t);

        // Residual and objective
        auto r = sub(C, P);
        double f = 0.5 * dot(r, r);

        // Gradient g = d/dt(0.5 * |r|^2) = r · C'
        double g = dot(r, Cp);

        // Hessian H = d^2/dt^2(0.5 * |r|^2) = |C'|^2 + r · C''
        double H = dot(Cp, Cp) + dot(r, Cpp);

        double gradNorm = std::abs(g);
        result.gradNorm = gradNorm;

        // Update best-so-far
        if (f < f_best)
        {
            f_best = f;
            t_best = t;
            C_best = C;
        }

        // Convergence by gradient norm
        if (gradNorm < tol)
        {
            bool finalOnBoundary =
                (t <= tol || t >= 1.0 - tol);

            result.status = finalOnBoundary
                                ? ClosestPointResult::Status::Boundary
                                : ClosestPointResult::Status::Success;
            break;
        }

        // Damped scalar Hessian: H + lambda
        double denom = H + lambda;

        // Hessian degeneracy / near-singularity
        if (std::abs(denom) < NEAR_ZERO)
        {
            // Gradient descent fallback
            double alpha = 0.1; // small step
            double delta_gd = -alpha * g;
            double t_gd = t + delta_gd;
            project_to_domain(t_gd);

            auto C_gd = evaluate(t_gd);
            auto r_gd = sub(C_gd, P);
            double f_gd = 0.5 * dot(r_gd, r_gd);

            if (f_gd < f)
            {
                // Accept GD step
                t = t_gd;

                // Update best-so-far
                if (f_gd < f_best)
                {
                    f_best = f_gd;
                    t_best = t_gd;
                    C_best = C_gd;
                }

                // Reset damping to encourage Newton-like behavior
                lambda = LAMBDA_SEED;
                continue;
            }

            // If GD didn't help, escalate damping
            lambda *= 10.0;
            if (lambda > LAMBDA_MAX)
            {
                result.status = ClosestPointResult::Status::Divergence;
                break;
            }

            continue;
        }

        // LM step: (H + lambda) * delta = -g
        double delta = -g / denom;
        double stepNorm = std::abs(delta);

        // Stagnation detection (tiny step + large damping)
        if (stepNorm < tol * 0.1 && lambda > LAMBDA_LARGE)
        {
            // Gradient descent fallback
            double alpha = 0.1;
            double delta_gd = -alpha * g;
            double t_gd = t + delta_gd;
            project_to_domain(t_gd);

            auto C_gd = evaluate(t_gd);
            auto r_gd = sub(C_gd, P);
            double f_gd = 0.5 * dot(r_gd, r_gd);

            if (f_gd < f)
            {
                // Accept GD step
                t = t_gd;

                // Update best-so-far
                if (f_gd < f_best)
                {
                    f_best = f_gd;
                    t_best = t_gd;
                    C_best = C_gd;
                }

                // Reset lambda to encourage Newton-like behavior
                lambda = LAMBDA_SEED;
                continue;
            }

            // GD didn't help -> declare stagnation
            result.status = ClosestPointResult::Status::Stagnation;
            break;
        }

        // Tentative LM update
        double t_new = t + delta;
        project_to_domain(t_new);

        auto C_new = evaluate(t_new);
        auto r_new = sub(C_new, P);
        double f_new = 0.5 * dot(r_new, r_new);

        // Update best-so-far
        if (f_new < f_best)
        {
            f_best = f_new;
            t_best = t_new;
            C_best = C_new;
        }

        // Accept or reject step
        if (f_new < f)
        {
            // Accept
            t = t_new;

            lambda *= 0.5;
            if (lambda < LAMBDA_MIN)
                lambda = LAMBDA_MIN;

            // Additional convergence check on step size
            if (stepNorm < tol)
            {
                bool finalOnBoundary =
                    (t <= tol || t >= 1.0 - tol);

                result.status = finalOnBoundary
                                    ? ClosestPointResult::Status::Boundary
                                    : ClosestPointResult::Status::Success;
                break;
            }
        }
        else
        {
            // Reject step, increase damping
            lambda *= 10.0;
            if (lambda > LAMBDA_MAX)
            {
                result.status = ClosestPointResult::Status::Divergence;
                break;
            }
        }
    }

    // Fill final result based on best-so-far
    result.iterations = iter;

    project_to_domain(t_best);
    auto C_final = evaluate(t_best);
    auto r_final = sub(C_final, P);
    double dist = norm(r_final);

    result.t() = t_best;
    result.point3D = C_final;
    result.distance = dist;
    result.bestObjective = f_best;

    result.onBoundary = (t_best <= tol || t_best >= 1.0 - tol);

    // If no status was set inside the loop, it means we exited by maxIters
    if (result.status != ClosestPointResult::Status::Success &&
        result.status != ClosestPointResult::Status::Boundary &&
        result.status != ClosestPointResult::Status::Stagnation &&
        result.status != ClosestPointResult::Status::Divergence)
    {
        result.status = ClosestPointResult::Status::MaxIterations;
    }

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