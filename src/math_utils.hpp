#pragma once
#include <array>
#include <cmath>
#include <stdexcept>

namespace MathUtils
{
    inline double dot(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    }

    inline std::array<double, 3> normalize(
        const std::array<double, 3> &v)
    {
        double len = std::sqrt(dot(v, v));
        if (len < 1e-14)
            return {0.0, 0.0, 0.0};

        return {v[0] / len, v[1] / len, v[2] / len};
    }

    inline std::array<double, 3> cross(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
    }

    inline std::array<double, 3> add(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return {a[0] + b[0],
                a[1] + b[1],
                a[2] + b[2]};
    }

    inline std::array<double, 3> sub(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b)
    {
        return {a[0] - b[0],
                a[1] - b[1],
                a[2] - b[2]};
    }

    inline std::array<double, 3> scale(
        const std::array<double, 3> &a,
        double s)
    {
        return {a[0] * s,
                a[1] * s,
                a[2] * s};
    }

    inline double norm(const std::array<double, 3> &a)
    {
        return std::sqrt(dot(a, a));
    }

    inline double WrapToInterval(double t, double t_min, double t_max)
    {
        const double L = t_max - t_min;

        if (L <= 0.0)
            throw std::runtime_error("WrapToInterval: invalid interval");

        double v = std::fmod(t - t_min, L);
        if (v < 0.0)
            v += L;
        return t_min + v;
    }
}