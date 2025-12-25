#include <array>
#include <cmath>

namespace MathUtils
{
    std::array<double, 3> normalize(
        const std::array<double, 3> &v);

    std::array<double, 3> cross(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b);

    double dot(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b);

    std::array<double, 3> add(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b);

    std::array<double, 3> sub(
        const std::array<double, 3> &a,
        const std::array<double, 3> &b);

    std::array<double, 3> scale(
        const std::array<double, 3> &a,
        double s);

    double norm(const std::array<double, 3> &a);
}