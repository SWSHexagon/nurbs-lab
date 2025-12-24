#include "nurbs_curve.hpp"
#include <cassert>
#include <iostream>

NURBSCurve::NURBSCurve(std::vector<std::array<double, 3>> control_points,
                       std::vector<double> weights,
                       std::vector<double> knots,
                       int degree)
    : m_ctrl(std::move(control_points)), m_w(std::move(weights)), m_knots(std::move(knots)),
      m_degree(degree), m_basis(m_degree, m_knots) // or whatever your BSplineBasis ctor is
{
    initialize_domain();
}

void NURBSCurve::initialize_domain()
{
    if (m_ctrl.empty())
        throw std::runtime_error("NURBSCurve: control point array is empty");

    if (m_ctrl.size() != m_w.size())
        throw std::runtime_error("NURBSCurve: weights must match control points");

    if (m_knots.size() != m_ctrl.size() + m_degree + 1)
        throw std::runtime_error("NURBSCurve: knot vector size mismatch");

    if (m_degree < 1 || m_degree > static_cast<int>(m_ctrl.size()) - 1)
        throw std::runtime_error("NURBSCurve: degree out of range");

    // Let the basis compute the domain (with padding)
    auto [umin, umax] = m_basis.getExtents();
    t_min = umin;
    t_max = umax;

    if (t_min >= t_max)
        throw std::runtime_error("NURBSCurve: invalid domain (t_min >= t_max)");
}

std::array<double, 3> NURBSCurve::evaluate(double t) const
{
    double tc = t;
    project_to_domain(tc);

    const int n = m_basis.numBasis(); // number of control points
    std::array<double, 3> S{0.0, 0.0, 0.0};
    double W = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double Ni = m_basis.evaluate(i, tc); // <-- YOUR API
        double Ni_wi = Ni * m_w[i];

        S = add(S, scale(m_ctrl[i], Ni_wi));
        W += Ni_wi;
    }

    if (W == 0.0)
        return {0.0, 0.0, 0.0};

    return scale(S, 1.0 / W);
}

std::array<double, 3> NURBSCurve::derivative(double t) const
{
    double tc = t;
    project_to_domain(tc);

    const int n = m_basis.numBasis();

    std::array<double, 3> S{0.0, 0.0, 0.0};
    std::array<double, 3> Sd{0.0, 0.0, 0.0};
    double W = 0.0;
    double Wd = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double Ni = m_basis.evaluate(i, tc);
        double dNi = m_basis.derivative(i, tc);

        double wi = m_w[i];

        S = add(S, scale(m_ctrl[i], Ni * wi));
        Sd = add(Sd, scale(m_ctrl[i], dNi * wi));

        W += Ni * wi;
        Wd += dNi * wi;
    }

    if (W == 0.0)
        return {0.0, 0.0, 0.0};

    double invW = 1.0 / W;
    double invW2 = invW * invW;

    auto num = sub(scale(Sd, W), scale(S, Wd));
    return scale(num, invW2);
}

std::array<double, 3> NURBSCurve::second_derivative(double t) const
{
    double tc = t;
    project_to_domain(tc);

    const int n = m_basis.numBasis();

    std::array<double, 3> S{0.0, 0.0, 0.0};
    std::array<double, 3> Sd{0.0, 0.0, 0.0};
    std::array<double, 3> Sdd{0.0, 0.0, 0.0};

    double W = 0.0;
    double Wd = 0.0;
    double Wdd = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double Ni = m_basis.evaluate(i, tc);
        double dNi = m_basis.derivative(i, tc);
        double ddNi = m_basis.second_derivative(i, tc);

        double wi = m_w[i];

        S = add(S, scale(m_ctrl[i], Ni * wi));
        Sd = add(Sd, scale(m_ctrl[i], dNi * wi));
        Sdd = add(Sdd, scale(m_ctrl[i], ddNi * wi));

        W += Ni * wi;
        Wd += dNi * wi;
        Wdd += ddNi * wi;
    }

    if (W == 0.0)
        return {0.0, 0.0, 0.0};

    double invW = 1.0 / W;
    double invW2 = invW * invW;
    double invW3 = invW2 * invW;

    // Numerator of C''
    auto term1 = scale(Sdd, W * W);
    auto term2 = scale(Sd, 2 * W * Wd);
    auto term3 = scale(S, W * Wdd);
    auto term4 = scale(S, 2 * Wd * Wd);

    auto num = add(sub(term1, term2), sub(term4, term3));

    return scale(num, invW3);
}

std::array<double, 3> NURBSCurve::tangent(double t, bool unitize) const
{
    auto d = derivative(t);
    return unitize ? normalize(d) : d;
}

std::array<double, 3> NURBSCurve::normal(double t) const
{
    // Same logic as CircleCurve, or a generic Frenet computation:
    auto T = tangent(t);
    auto dT = derivative(t); // or use second_derivative if you prefer

    auto proj = scale(T, dot(T, dT));
    auto Np = sub(dT, proj);
    return normalize(Np);
}

void NURBSCurve::DumpInfo() const
{
    std::cout << "NURBSCurve:\n";
    std::cout << "  Degree: " << m_degree << "\n";
    std::cout << "  Control points: " << m_ctrl.size() << "\n";
    std::cout << "  Domain: [" << t_min << ", " << t_max << "]\n";
}
