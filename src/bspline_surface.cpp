#include "bspline_surface.hpp"
#include "math_utils.hpp"

using namespace MathUtils;

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

std::array<double, 3> BSplineSurface::derivative_u(double u, double v) const
{
    std::array<double, 3> Su = {0, 0, 0};

    for (int i = 0; i < u_basis_.numBasis(); ++i)
    {
        double Ni_u = u_basis_.derivative(i, u);

        for (int j = 0; j < v_basis_.numBasis(); ++j)
        {
            double Mj = v_basis_.evaluate(j, v);

            const auto &P = ctrl_[i][j];

            Su[0] += Ni_u * Mj * P[0];
            Su[1] += Ni_u * Mj * P[1];
            Su[2] += Ni_u * Mj * P[2];
        }
    }

    return Su;
}

std::array<double, 3> BSplineSurface::derivative_v(double u, double v) const
{
    std::array<double, 3> Sv = {0, 0, 0};

    for (int i = 0; i < u_basis_.numBasis(); ++i)
    {
        double Ni = u_basis_.evaluate(i, u);

        for (int j = 0; j < v_basis_.numBasis(); ++j)
        {
            double Mj_v = v_basis_.derivative(j, v);

            const auto &P = ctrl_[i][j];

            Sv[0] += Ni * Mj_v * P[0];
            Sv[1] += Ni * Mj_v * P[1];
            Sv[2] += Ni * Mj_v * P[2];
        }
    }

    return Sv;
}

std::array<double, 3> BSplineSurface::second_derivative_uu(double u, double v) const
{
    std::array<double, 3> Suu = {0, 0, 0};

    for (int i = 0; i < u_basis_.numBasis(); ++i)
    {
        double Ni_uu = u_basis_.second_derivative(i, u);

        for (int j = 0; j < v_basis_.numBasis(); ++j)

        {
            double Mj = v_basis_.evaluate(j, v);
            const auto &P = ctrl_[i][j];

            Suu[0] += Ni_uu * Mj * P[0];
            Suu[1] += Ni_uu * Mj * P[1];
            Suu[2] += Ni_uu * Mj * P[2];
        }
    }

    return Suu;
}

std::array<double, 3> BSplineSurface::second_derivative_vv(double u, double v) const
{
    std::array<double, 3> Svv = {0, 0, 0};

    for (int i = 0; i < u_basis_.numBasis(); ++i)
    {
        double Ni = u_basis_.evaluate(i, u);

        for (int j = 0; j < v_basis_.numBasis(); ++j)
        {
            double Mj_vv = v_basis_.second_derivative(j, v);
            const auto &P = ctrl_[i][j];

            Svv[0] += Ni * Mj_vv * P[0];
            Svv[1] += Ni * Mj_vv * P[1];
            Svv[2] += Ni * Mj_vv * P[2];
        }
    }

    return Svv;
}

std::array<double, 3> BSplineSurface::second_derivative_uv(double u, double v) const
{
    std::array<double, 3> Suv = {0, 0, 0};

    for (int i = 0; i < u_basis_.numBasis(); ++i)
    {
        double Ni_u = u_basis_.derivative(i, u);

        for (int j = 0; j < v_basis_.numBasis(); ++j)
        {
            double Mj_v = v_basis_.derivative(j, v);
            const auto &P = ctrl_[i][j];

            Suv[0] += Ni_u * Mj_v * P[0];
            Suv[1] += Ni_u * Mj_v * P[1];
            Suv[2] += Ni_u * Mj_v * P[2];
        }
    }

    return Suv;
}

std::array<double, 3> BSplineSurface::normal(double u, double v) const
{
    auto Su = derivative_u(u, v);
    auto Sv = derivative_v(u, v);

    std::array<double, 3> n;

    n[0] = Su[1] * Sv[2] - Su[2] * Sv[1];
    n[1] = Su[2] * Sv[0] - Su[0] * Sv[2];
    n[2] = Su[0] * Sv[1] - Su[1] * Sv[0];

    return (normalize(n));
}

inline void BSplineSurface::project_to_domain(double &u, double &v)
{
    // If already inside the domain, nothing to do
    if (u >= 0.0 && u <= 1.0 &&
        v >= 0.0 && v <= 1.0)
    {
        return;
    }

    // Case 1: Outside on the left or right → project to vertical edges
    if (u < 0.0)
        u = 0.0;
    else if (u > 1.0)
        u = 1.0;

    // Case 2: Outside on the bottom or top → project to horizontal edges
    if (v < 0.0)
        v = 0.0;
    else if (v > 1.0)
        v = 1.0;

    // At this point, (u,v) is guaranteed to be in the closed domain.
    // The logic above naturally handles:
    //   - corner projection
    //   - edge projection
    //   - interior no-op
}

BSplineSurface::ClosestPointResult BSplineSurface::closest_point_LM(
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
        auto Su = derivative_u(u, v);
        auto Sv = derivative_v(u, v);

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

    result.onBoundary = (u_best <= tol || u_best >= 1.0 - tol ||
                         v_best <= tol || v_best >= 1.0 - tol);

    return result;
}

BSplineSurface::SurfaceDifferential BSplineSurface::differential(double u, double v) const
{
    SurfaceDifferential d;

    d.Su = derivative_u(u, v);
    d.Sv = derivative_v(u, v);
    d.Suu = second_derivative_uu(u, v);
    d.Suv = second_derivative_uv(u, v);
    d.Svv = second_derivative_vv(u, v);

    auto n_raw = cross(d.Su, d.Sv);
    d.N = normalize(n_raw);

    d.E = dot(d.Su, d.Su);
    d.F = dot(d.Su, d.Sv);
    d.G = dot(d.Sv, d.Sv);

    d.e = dot(d.Suu, d.N);
    d.f = dot(d.Suv, d.N);
    d.g = dot(d.Svv, d.N);

    return d;
}

BSplineSurface::Curvature BSplineSurface::curvature(double u, double v) const
{
    Curvature c{0, 0, 0, 0, false};

    auto d = differential(u, v);

    double denom = d.E * d.G - d.F * d.F;

    if (std::abs(denom) < 1e-14)
        return c; // degenerate metric (e.g. collapsed patch)

    c.K = (d.e * d.g - d.f * d.f) / denom;
    c.H = (d.E * d.g - 2.0 * d.F * d.f + d.G * d.e) / (2.0 * denom);

    double disc = c.H * c.H - c.K;

    if (disc < 0.0)
        disc = 0.0; // clamp small negative from numerical noise

    double root = std::sqrt(disc);

    c.k1 = c.H + root;
    c.k2 = c.H - root;
    c.valid = true;

    return c;
}

void BSplineSurface::DumpInfo() const
{
    u_basis_.DumpInfo("U Basis Knots:");
    v_basis_.DumpInfo("V Basis Knots:");
}
