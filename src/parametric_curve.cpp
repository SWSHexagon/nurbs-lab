#include "parametric_curve.hpp"
#include "math_utils.hpp"
#include <cmath>
#include <algorithm>

using namespace MathUtils;

void ParametricCurve::project_to_domain(double &t) const
{
    if (m_is_periodic)
    {
        t = WrapToInterval(t, t_min, t_max);
        return;
    }

    t = std::clamp(t, t_min, t_max);
}

ClosestPointResult ParametricCurve::closest_point_LM(
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
            if (!m_is_periodic)
            {
                bool finalOnBoundary = (t < t_min + tol || t > t_max - tol);

                result.status = finalOnBoundary
                                    ? ClosestPointResult::Status::Boundary
                                    : ClosestPointResult::Status::Success;
            }
            else
                result.status = ClosestPointResult::Status::Success;

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
                bool finalOnBoundary = (t < t_min + tol || t > t_max - tol);

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
    auto Cp_final = derivative(t_best);
    auto r_final = sub(C_final, P);
    double dist = norm(r_final);

    // Set gradient at t_best
    double g_final = dot(r_final, Cp_final);
    double gradNorm_final = std::abs(abs(g_final));
    result.gradNorm = gradNorm_final;

    result.t() = t_best;
    result.point3D = C_final;
    result.distance = dist;
    result.bestObjective = f_best;

    result.onBoundary = (t_best < t_min + tol || t_best > t_max - tol);

    // Stagnation at a point with gradient essentially zero is success
    if (result.status == ClosestPointResult::Status::Stagnation &&
        gradNorm_final < tol * 10)
    {
        if (!m_is_periodic)
        {
            result.status = result.onBoundary
                                ? ClosestPointResult::Status::Boundary
                                : ClosestPointResult::Status::Success;
        }
        else
            result.status = ClosestPointResult::Status::Success;
    }

    return result;
}
