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
    // Initialize surface extents
    initialize_domain();

    // Build the boundary curves
    build_boundary_curves();
}

void BSplineSurface::initialize_domain()
{
    // Initialize surface extents
    auto [uMin, uMax] = u_basis_.getExtents();
    auto [vMin, vMax] = v_basis_.getExtents();

    u_min = uMin;
    u_max = uMax;
    v_min = vMin;
    v_max = vMax;
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

std::array<double, 3> BSplineSurface::normal(double u, double v, bool unitize /*= true*/) const
{
    auto Su = derivative_u(u, v);
    auto Sv = derivative_v(u, v);

    std::array<double, 3> n;

    n[0] = Su[1] * Sv[2] - Su[2] * Sv[1];
    n[1] = Su[2] * Sv[0] - Su[0] * Sv[2];
    n[2] = Su[0] * Sv[1] - Su[1] * Sv[0];

    if (!unitize)
        return (n);
    else
        return (normalize(n));
}

void BSplineSurface::build_boundary_curves()
{
    const int nU = ctrl_.size();    // num_ctrl_u();
    const int nV = ctrl_[0].size(); // num_ctrl_v();
    const int degree_u = u_basis_.degree();
    const int degree_v = v_basis_.degree();

    // 1) u = u_min edge: varies in v, so use v-basis
    {
        BSplineBasis basis_v(degree_v, v_basis_.knots());
        std::vector<std::array<double, 3>> cps(nV);

        int i = 0; // left column
        for (int j = 0; j < nV; ++j)
            cps[j] = ctrl_[i][j];

        boundary_u_min_ = BSplineCurve(basis_v, std::move(cps));
    }

    // 2) u = u_max edge: varies in v, so use v-basis
    {
        BSplineBasis basis_v(degree_v, v_basis_.knots());
        std::vector<std::array<double, 3>> cps(nV);

        int i = nU - 1; // right column
        for (int j = 0; j < nV; ++j)
            cps[j] = ctrl_[i][j];

        boundary_u_max_ = BSplineCurve(basis_v, std::move(cps));
    }

    // 3) v = v_min edge: varies in u, so use u-basis
    {
        BSplineBasis basis_u(degree_u, u_basis_.knots());
        std::vector<std::array<double, 3>> cps(nU);

        int j = 0; // bottom row
        for (int i = 0; i < nU; ++i)
            cps[i] = ctrl_[i][j];

        boundary_v_min_ = BSplineCurve(basis_u, std::move(cps));
    }

    // 4) v = v_max edge: varies in u, so use u-basis
    {
        BSplineBasis basis_u(degree_u, u_basis_.knots());
        std::vector<std::array<double, 3>> cps(nU);

        int j = nV - 1; // top row
        for (int i = 0; i < nU; ++i)
            cps[i] = ctrl_[i][j];

        boundary_v_max_ = BSplineCurve(basis_u, std::move(cps));
    }
}

std::array<double, 4> BSplineSurface::domain() const
{
    return std::array<double, 4>{u_min, u_max, v_min, v_max};
}

SurfaceDifferential BSplineSurface::differential(double u, double v) const
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

    d.metricDet = d.E * d.G - d.F * d.F;
    d.metricDegenerate = (d.metricDet < NEAR_ZERO);

    return d;
}

SurfaceCurvature BSplineSurface::curvature(double u, double v) const
{
    SurfaceCurvature c{0, 0, 0, 0, false};

    auto d = differential(u, v);

    if (d.metricDegenerate)
        return (c); // degenerate metric (e.g. collapsed patch)

    double denom = d.metricDet;
    c.K = (d.e * d.g - d.f * d.f) / denom;
    c.H = (d.E * d.g - 2.0 * d.F * d.f + d.G * d.e) / (2.0 * denom);

    double disc = c.H * c.H - c.K;

    if (disc < 0.0)
        disc = 0.0; // clamp small negative from numerical noise

    double root = std::sqrt(disc);

    c.k1 = c.H + root;
    c.k2 = c.H - root;
    c.valid = true;

    // --- Principal directions in UV ---
    auto uv1 = principal_direction_uv(d, c.k1);
    auto uv2 = principal_direction_uv(d, c.k2);

    c.uv1 = uv1;
    c.uv2 = uv2;

    // --- Map to 3D ---
    c.dir1 = combine_uv_to_3d(d, uv1);
    c.dir2 = combine_uv_to_3d(d, uv2);

    return (c);
}

std::array<double, 2> BSplineSurface::principal_direction_uv(
    const SurfaceDifferential &d,
    double k)
{
    double a = d.e - k * d.E;
    double b = d.f - k * d.F;
    double c = d.g - k * d.G;

    std::array<double, 2> uv{0.0, 0.0};

    // Choose a stable null-space vector for [a b; b c]
    if (std::abs(a) + std::abs(b) > std::abs(b) + std::abs(c))
    {
        // Use (-b, a)
        uv[0] = -b;
        uv[1] = a;
    }
    else
    {
        // Use (-c, b)
        uv[0] = -c;
        uv[1] = b;
    }

    double len = std::sqrt(uv[0] * uv[0] + uv[1] * uv[1]);
    if (len < NEAR_ZERO)
    {
        // Fallback: arbitrary direction
        uv[0] = 1.0;
        uv[1] = 0.0;
    }
    else
    {
        uv[0] /= len;
        uv[1] /= len;
    }

    return (uv);
}

std::array<double, 3> BSplineSurface::combine_uv_to_3d(
    const SurfaceDifferential &d,
    const std::array<double, 2> &uv)
{
    std::array<double, 3> r{
        uv[0] * d.Su[0] + uv[1] * d.Sv[0],
        uv[0] * d.Su[1] + uv[1] * d.Sv[1],
        uv[0] * d.Su[2] + uv[1] * d.Sv[2]};

    double len = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    if (len < 1e-14)
        return {0.0, 0.0, 0.0};

    r[0] /= len;
    r[1] /= len;
    r[2] /= len;

    return (r);
}

SurfaceFrame BSplineSurface::frame(double u, double v) const
{
    SurfaceFrame fr;

    // Position
    fr.P = evaluate(u, v);

    // Differential geometry
    auto d = differential(u, v);

    fr.Su = d.Su;
    fr.Sv = d.Sv;
    fr.N = d.N;

    fr.E = d.E;
    fr.F = d.F;
    fr.G = d.G;

    fr.e = d.e;
    fr.f = d.f;
    fr.g = d.g;

    fr.metricDet = d.metricDet;
    fr.metricDegenerate = d.metricDegenerate;

    // Curvature
    auto c = curvature(u, v);
    fr.k1 = c.k1;
    fr.k2 = c.k2;
    fr.valid = c.valid;

    fr.dir1 = c.dir1;
    fr.dir2 = c.dir2;

    fr.uv1 = c.uv1;
    fr.uv2 = c.uv2;

    return (fr);
}

void BSplineSurface::DumpInfo() const
{
    u_basis_.DumpInfo("U Basis Knots:");
    v_basis_.DumpInfo("V Basis Knots:");
}
