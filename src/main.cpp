#include "bspline_surface_builder.hpp"
#include <fstream>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
    std::cout << "Running nurbs_demo...\n";

    // ============================================================
    // 1. BASIS FUNCTION CSV OUTPUT (degree 2, 4 basis functions)
    // ============================================================
    BSplineBasis basis(2, {0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1}); // 8 knots → 4 basis functions
    basis.GeneratePlot("plots/data/basis.csv");

    std::array<double, 3> p0 = {2, 3, 0};
    std::array<double, 3> u_dir = {1, 0, 0};
    std::array<double, 3> v_dir = {0, 1, 0};

    auto surf = SurfaceBuilder::Plane(1, 4, p0, u_dir, v_dir);

    // ============================================================
    // 3. SAMPLE SURFACE AND EXPORT
    // ============================================================
    std::ofstream sout("plots/data/surface.csv");
    std::ofstream tout("plots/data/tsurface.xyz"); // XYZ for xyzviewer.com

    if (!sout || !tout)
    {
        std::cerr << "Failed to open output files\n";
        return 1;
    }

    int ptCount = 0;

    for (int i = 0; i <= 20; ++i)
    {
        double u = double(i) / 20.0;

        if (u > 1.0)
            u = 1.0;

        for (int j = 0; j <= 20; ++j)
        {
            double v = double(j) / 20.0;

            if (v > 1.0)
                v = 1.0;

            auto p = surf.evaluate(u, v);

            // CSV for plotting
            sout << u << "," << v << "," << p[0] << "," << p[1] << "," << p[2] << "\n";

            // XYZ for point cloud viewer
            tout << p[0] << " " << p[1] << " " << p[2] << "\n";

            ptCount++;
        }
    }

    std::cout << "Demo complete.  Wrote " << ptCount << " points.\n";
    return 0;
}
