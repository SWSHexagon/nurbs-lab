#include "bspline_surface.hpp"
#include <iostream>
#include <cmath>

void test_sphere_curvature(const BSplineSurface &surf, double R)
{
    std::cout << "=== Sphere Curvature Test ===\n";

    // Sample a grid of UV points
    for (int iu = 0; iu <= 5; ++iu)
    {
        for (int iv = 0; iv <= 5; ++iv)
        {
            double u = iu / 5.0;
            double v = iv / 5.0;

            auto c = surf.curvature(u, v);

            if (!c.valid)
            {
                std::cout << "Invalid differential at (" << u << "," << v << ")\n";
                continue;
            }

            double k_expected = 1.0 / R;

            std::cout << "u=" << u << " v=" << v
                      << "  k1=" << c.k1
                      << "  k2=" << c.k2
                      << "  expected=" << k_expected << "\n";

            // Check closeness
            if (std::abs(c.k1 - k_expected) > 1e-3 ||
                std::abs(c.k2 - k_expected) > 1e-3)
            {
                std::cout << "  **FAIL** curvature mismatch\n";
            }
        }
    }
}

int main()
{
    std::cout << "=== Testing sphere curvature ===\n";
    return 0;
}