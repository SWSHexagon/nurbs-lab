#include "circle_curve.hpp"
#include <cmath>
#include <numbers>
#include <iostream>

CircleCurve::CircleCurve(
    const std::array<double, 3> &center,
    const std::array<double, 3> &normal,
    double radius)
    : C(center), R(radius)
{
    N = normalize(normal);

    // Build orthonormal basis (U, V) in the plane of the circle
    std::array<double, 3> tmp = {1, 0, 0};
    if (std::abs(dot(tmp, N)) > 0.9)
        tmp = {0, 1, 0};

    U = normalize(sub(tmp, scale(N, dot(tmp, N))));
    V = normalize({N[1] * U[2] - N[2] * U[1],
                   N[2] * U[0] - N[0] * U[2],
                   N[0] * U[1] - N[1] * U[0]});

    initialize_domain();
}

void CircleCurve::initialize_domain()
{
    t_min = 0.0;
    t_max = 2.0 * std::numbers::pi;
}

std::array<double, 3> CircleCurve::evaluate(double t) const
{
    return add(C, add(scale(U, R * std::cos(t)),
                      scale(V, R * std::sin(t))));
}

std::array<double, 3> CircleCurve::derivative(double t) const
{
    return add(scale(U, -R * std::sin(t)),
               scale(V, R * std::cos(t)));
}

std::array<double, 3> CircleCurve::second_derivative(double t) const
{
    return add(scale(U, -R * std::cos(t)),
               scale(V, -R * std::sin(t)));
}

std::array<double, 3> CircleCurve::tangent(double t, bool unitize) const
{
    auto d = derivative(t);
    return unitize ? normalize(d) : d;
}

std::array<double, 3> CircleCurve::normal(double t) const
{
    // Principal normal = derivative of unit tangent
    auto T = tangent(t);
    auto dT = derivative(t);
    auto proj = scale(T, dot(T, dT));
    return normalize(sub(dT, proj));
}

void CircleCurve::DumpInfo() const
{
    std::cout << "CircleCurve:\n";
    std::cout << "  Center: (" << C[0] << ", " << C[1] << ", " << C[2] << ")\n";
    std::cout << "  Radius: " << R << "\n";
    std::cout << "  Domain: [" << t_min << ", " << t_max << "]\n";
}