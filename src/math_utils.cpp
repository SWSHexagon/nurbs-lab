#include "math_utils.hpp"

namespace MathUtils
{
    std::array<double, 3> normalize(
        const std::array<double, 3> &v)
    {
        double len = std::sqrt(dot(v, v));
        if (len < 1e-14)
            return {0.0, 0.0, 0.0};

        return {v[0] / len, v[1] / len, v[2] / len};
    }

    std::array<double, 3> cross(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
    }

    double dot(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    }

    std::array<double, 3> add(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return {a[0] + b[0],
                a[1] + b[1],
                a[2] + b[2]};
    }

    std::array<double, 3> sub(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return {a[0] - b[0],
                a[1] - b[1],
                a[2] - b[2]};
    }

    std::array<double, 3> scale(
        const std::array<double, 3> &a,
        double s)
    {
        return {a[0] * s,
                a[1] * s,
                a[2] * s};
    }

    double norm(const std::array<double, 3> &a)
    {
        return std::sqrt(dot(a, a));
    }
} // namespace MathUtils