#include "parametric_surface.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <boost/asio.hpp>
#include <boost/thread.hpp>

void ParametricSurface::project_to_domain(double &u, double &v) const
{
    u = std::clamp(u, u_min, u_max);
    v = std::clamp(v, v_min, v_max);
}

ClosestPointResult ParametricSurface::closest_point_LM(
    const std::array<double, 3> &point,
    double u0,
    double v0,
    int maxIter,
    double tol) const
{
    ClosestPointResult result;

    double u = u0;
    double v = v0;

    double u_best = u;
    double v_best = v;

    auto S0 = evaluate(u, v);
    std::array<double, 3> F0 = {
        S0[0] - point[0],
        S0[1] - point[1],
        S0[2] - point[2]};
    double f_best = 0.5 * (F0[0] * F0[0] + F0[1] * F0[1] + F0[2] * F0[2]);

    double lambda = LAMBDA_SEED;
    int iter = 0;

    for (iter = 0; iter < maxIter; ++iter)
    {
        // Evaluate surface and derivatives
        auto S = evaluate(u, v);
        auto [Su, Sv] = derivatives(u, v);

        // Residual F = S - P
        std::array<double, 3> F = {
            S[0] - point[0],
            S[1] - point[1],
            S[2] - point[2]};

        // Gradient of 0.5||F||^2
        double Fu = Su[0] * F[0] + Su[1] * F[1] + Su[2] * F[2];
        double Fv = Sv[0] * F[0] + Sv[1] * F[1] + Sv[2] * F[2];

        double gradNorm = std::sqrt(Fu * Fu + Fv * Fv);
        result.gradNorm = gradNorm;

        // Check for convergence at an extrema or boundary
        if (gradNorm < tol)
        {
            double f_here = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);

            if (f_here < f_best)
            {
                f_best = f_here;
                u_best = u;
                v_best = v;
            }

            bool finalOnBoundary =
                (u <= tol || u >= 1.0 - tol ||
                 v <= tol || v >= 1.0 - tol);

            result.status = finalOnBoundary
                                ? ClosestPointResult::Status::Boundary
                                : ClosestPointResult::Status::Success;

            break;
        }

        // Gauss–Newton Hessian approximation
        double Huu = Su[0] * Su[0] + Su[1] * Su[1] + Su[2] * Su[2];
        double Hvv = Sv[0] * Sv[0] + Sv[1] * Sv[1] + Sv[2] * Sv[2];
        double Huv = Su[0] * Sv[0] + Su[1] * Sv[1] + Su[2] * Sv[2];

        // Damped Hessian
        double a = Huu + lambda;
        double b = Huv;
        double c = Huv;
        double d = Hvv + lambda;

        // Determinant check
        double det = a * d - b * c;
        double scale = std::max({std::abs(a), std::abs(b), std::abs(c), std::abs(d), 1.0});

        if (std::abs(det) < NEAR_ZERO * scale * scale)
        {
            // Gradient descent fallback
            double alpha = 0.1; // small step
            double du_gd = -alpha * Fu;
            double dv_gd = -alpha * Fv;

            double u_gd = u + du_gd;
            double v_gd = v + dv_gd;
            project_to_domain(u_gd, v_gd);

            auto Sgd = evaluate(u_gd, v_gd);
            std::array<double, 3> Fgd = {
                Sgd[0] - point[0],
                Sgd[1] - point[1],
                Sgd[2] - point[2]};

            double f_old = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
            double f_gd = 0.5 * (Fgd[0] * Fgd[0] + Fgd[1] * Fgd[1] + Fgd[2] * Fgd[2]);

            if (f_gd < f_old)
            {
                // Accept fallback step
                u = u_gd;
                v = v_gd;

                // Update best-so-far as well
                if (f_gd < f_best)
                {
                    f_best = f_gd;
                    u_best = u_gd;
                    v_best = v_gd;
                }

                // Reset damping to encourage Newton-like behavior
                lambda = LAMBDA_SEED;
                continue;
            }

            // If GD didn't help, escalate damping as before
            lambda *= 10.0;
            if (lambda > LAMBDA_MAX)
            {
                result.status = ClosestPointResult::Status::Divergence;
                break;
            }
            continue;
        }

        // Solve for LM step
        double du = (-Fu * d + Fv * b) / det;
        double dv = (-Fv * a + Fu * c) / det;

        double stepNorm = std::sqrt(du * du + dv * dv);

        // Stagnation detection
        if (stepNorm < tol * 0.1 && lambda > LAMBDA_LARGE)
        {
            // 1. Try gradient descent fallback
            double alpha = 0.1;
            double du_gd = -alpha * Fu;
            double dv_gd = -alpha * Fv;

            double u_gd = u + du_gd;
            double v_gd = v + dv_gd;
            project_to_domain(u_gd, v_gd);

            auto Sgd = evaluate(u_gd, v_gd);
            std::array<double, 3> Fgd = {
                Sgd[0] - point[0],
                Sgd[1] - point[1],
                Sgd[2] - point[2]};

            double f_old = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
            double f_gd = 0.5 * (Fgd[0] * Fgd[0] + Fgd[1] * Fgd[1] + Fgd[2] * Fgd[2]);

            if (f_gd < f_old)
            {
                // Accept gradient descent step
                u = u_gd;
                v = v_gd;

                // Update best-so-far as well
                if (f_gd < f_best)
                {
                    f_best = f_gd;
                    u_best = u_gd;
                    v_best = v_gd;
                }

                // Reset lambda to encourage Newton-like behavior
                lambda = LAMBDA_SEED;
                continue;
            }

            // 2. If GD didn't help, declare stagnation
            result.status = ClosestPointResult::Status::Stagnation;
            break;
        }

        // Tentative update
        double u_new = u + du;
        double v_new = v + dv;

        // Project to domain
        project_to_domain(u_new, v_new);

        // Evaluate new residual
        auto Snew = evaluate(u_new, v_new);
        std::array<double, 3> Fnew = {
            Snew[0] - point[0],
            Snew[1] - point[1],
            Snew[2] - point[2]};

        double f_old = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
        double f_new = 0.5 * (Fnew[0] * Fnew[0] + Fnew[1] * Fnew[1] + Fnew[2] * Fnew[2]);

        if (f_new < f_best)
        {
            f_best = f_new;
            u_best = u_new;
            v_best = v_new;
        }

        // Accept or reject
        if (f_new < f_old)
        {
            u = u_new;
            v = v_new;

            lambda *= 0.5;
            if (lambda < LAMBDA_MIN)
                lambda = LAMBDA_MIN;

            if (stepNorm < tol)
            {
                bool finalOnBoundary =
                    (u <= tol || u >= 1.0 - tol ||
                     v <= tol || v >= 1.0 - tol);

                result.status = finalOnBoundary
                                    ? ClosestPointResult::Status::Boundary
                                    : ClosestPointResult::Status::Success;

                break;
            }
        }
        else
        {
            lambda *= 10.0;
            if (lambda > LAMBDA_MAX)
            {
                result.status = ClosestPointResult::Status::Divergence;
                break;
            }
        }
    }

    // Fill final result
    result.iterations = iter;
    result.u = u_best;
    result.v = v_best;

    auto Sfinal = evaluate(u_best, v_best);
    result.point3D = Sfinal;
    double dx = Sfinal[0] - point[0];
    double dy = Sfinal[1] - point[1];
    double dz = Sfinal[2] - point[2];
    result.distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    result.bestObjective = f_best;

    result.onBoundary = (u_best <= tol || u_best >= 1.0 - tol ||
                         v_best <= tol || v_best >= 1.0 - tol);

    return result;
}

ClosestPointResult ParametricSurface::closest_point_with_boundary_fallback(
    const std::array<double, 3> &P,
    double u0,
    double v0,
    int maxIter,
    double tol,
    int maxCurveIter,
    double curveTol) const
{
    // 1. Run the surface LM solver
    auto surfCP = closest_point_LM(P, u0, v0, maxIter, tol);

    // Start with surface result as current best
    ClosestPointResult best = surfCP;

    // 2. If not on boundary, nothing else to do
    if (!surfCP.onBoundary)
        return best;

    // Collect only the boundaries LM actually touched
    std::vector<const ParametricCurve *> candidates;

    if (surfCP.u <= tol)
        candidates.push_back(&boundary_u_min());
    if (surfCP.u >= 1.0 - tol)
        candidates.push_back(&boundary_u_max());
    if (surfCP.v <= tol)
        candidates.push_back(&boundary_v_min());
    if (surfCP.v >= 1.0 - tol)
        candidates.push_back(&boundary_v_max());

    // Example: try each with t0 = 0.5
    for (const ParametricCurve *c : candidates)
    {
        auto cp_curve = c->closest_point_LM(P, 0.5);
        if (cp_curve.distance < best.distance)
        {
            // translate cp_curve into a surface-style ClosestPointResult
            // using a (u,v) mapping associated with each boundary (next step)
            best.point3D = cp_curve.point3D;
            best.distance = cp_curve.distance;
            best.onBoundary = true;
            best.status = ClosestPointResult::Status::Boundary;

            if (c == &boundary_u_min())
            {
                best.u = 0.0;
                best.v = cp_curve.t();
            }
            else if (c == &boundary_u_max())
            {
                best.u = 1.0;
                best.v = cp_curve.t();
            }
            else if (c == &boundary_v_min())
            {
                best.u = cp_curve.t();
                best.v = 0.0;
            }
            else if (c == &boundary_v_max())
            {
                best.u = cp_curve.t();
                best.v = 1.0;
            }
        }
    }

    return best;
}

ClosestPointResult ParametricSurface::closest_point_global(
    const std::array<double, 3> &P,
    int gridResolution,
    double u0,
    double v0,
    int maxIter,
    double tol,
    int maxCurveIter,
    double curveTol) const
{
    boost::asio::thread_pool pool(boost::thread::hardware_concurrency());

    std::mutex bestMutex;
    ClosestPointResult best;
    best.distance = std::numeric_limits<double>::infinity();

    for (int i = 0; i < gridResolution; ++i)
    {
        double u_seed = double(i) / (gridResolution - 1);

        for (int j = 0; j < gridResolution; ++j)
        {
            double v_seed = double(j) / (gridResolution - 1);

            boost::asio::post(pool, [&, u_seed, v_seed]()
                              {
                auto result = closest_point_with_boundary_fallback(
                    P, u_seed, v_seed, maxIter, tol, maxCurveIter, curveTol);

                std::lock_guard<std::mutex> lock(bestMutex);
                if (result.distance < best.distance)
                    best = result; });
        }
    }

    pool.join();

    return best;
}
