#include "bspline_basis.hpp"
#include "bspline_surface.hpp"
#include "test_closest_point.hpp"
#include <fstream>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

void check_basis_consistency(const BSplineBasis &B, int numCtrl)
{
    const int p = B.degree();
    const auto &U = B.knots();
    const int m = static_cast<int>(U.size()) - 1; // last knot index
    const int nBasis = B.numBasis();              // expected basis count
    const int expectedBasis = static_cast<int>(U.size()) - p - 1;

    std::cout << "\n=== B-Spline Basis Consistency Check ===\n";

    // 1. Degree sanity
    if (p < 0)
        std::cout << "ERROR: degree < 0\n";

    // 2. Knot vector length sanity
    if (U.size() < static_cast<size_t>(p + 2))
        std::cout << "ERROR: knot vector too short for degree\n";

    // 3. Basis count sanity
    std::cout << "Basis functions: " << nBasis
              << " (expected " << expectedBasis << ")\n";

    if (nBasis != expectedBasis)
        std::cout << "ERROR: numBasis() mismatch\n";

    // 4. Control point count consistency
    if (numCtrl != nBasis)
        std::cout << "ERROR: control point count (" << numCtrl
                  << ") does not match basis count (" << nBasis << ")\n";

    // 5. Parameter domain
    double umin = U[p];
    double umax = U[m - p];
    std::cout << "Parameter domain: [" << umin << ", " << umax << "]\n";

    // 6. Endpoint behavior
    std::cout << "\nEndpoint basis values:\n";
    std::cout << "At u = " << umin << ":\n";
    for (int i = 0; i < nBasis; ++i)
        std::cout << "  N" << i << " = " << B.evaluate(i, umin) << "\n";

    std::cout << "At u = " << umax << ":\n";
    for (int i = 0; i < nBasis; ++i)
        std::cout << "  N" << i << " = " << B.evaluate(i, umax) << "\n";

    // 7. Partition of unity test
    std::cout << "\nPartition of unity test:\n";
    for (double u = umin; u <= umax; u += (umax - umin) / 10.0)
    {
        double sum = 0.0;
        for (int i = 0; i < nBasis; ++i)
            sum += B.evaluate(i, u);

        std::cout << "  u=" << std::setw(5) << u
                  << "  sum(Ni)=" << sum;

        if (std::abs(sum - 1.0) > 1e-6)
            std::cout << "  <-- ERROR";

        std::cout << "\n";
    }

    std::cout << "=== End Consistency Check ===\n\n";
}

int main()
{
    std::cout << "Running nurbs_demo...\n";

    // ============================================================
    // 1. BASIS FUNCTION CSV OUTPUT (degree 2, 4 basis functions)
    // ============================================================
    BSplineBasis basis(2, {0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1}); // 8 knots → 4 basis functions

    check_basis_consistency(basis, 5);

    std::ofstream bout("plots/data/basis.csv");
    if (!bout)
    {
        std::cerr << "Failed to open basis.csv\n";
        return 1;
    }

    for (int i = 0; i <= 100; i++)
    {
        double u = 0.01 * double(i);

        if (u > 1.0)
            u = 1.0;

        bout << u;
        for (int i = 0; i < basis.numBasis(); ++i)
            bout << "," << basis.evaluate(i, u);
        bout << "\n";
    }

    std::cout << "Wrote plots/data/basis.csv\n";

    // ============================================================
    // 2. B-SPLINE SURFACE SETUP (4x4 control net)
    // ============================================================
    BSplineBasis ub(2, {0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1});
    BSplineBasis vb(2, {0, 0, 0, 1.0 / 3.0, 2.0 / 3.0, 1, 1, 1});

    std::cout << "Basis check:\n";
    for (int i = 0; i < ub.numBasis(); ++i)
        std::cout << "N" << i << "(0) = " << ub.evaluate(i, 0.0) << "\n";
    for (int i = 0; i < ub.numBasis(); ++i)
        std::cout << "N" << i << "(1) = " << ub.evaluate(i, 1.0) << "\n";

    // 4×4 tilted curved control net
    std::vector<std::vector<std::array<double, 3>>> ctrl(5, std::vector<std::array<double, 3>>(5));

    const double ax = 0;   // tilt in x
    const double ay = 0;   // tilt in y
    const double k = 0.25; // curvature strength

    std::ofstream ctrlnet("plots/data/ctrl_net.xyz");
    if (!ctrlnet)
    {
        std::cerr << "Failed to open ctrl_net\n";
        return 1;
    }

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            double x = double(i) - 2.0;
            double y = double(j) - 2.0;
            double z = k * (x * x + y * y) + ax * x + ay * y;
            ctrl[i][j] = {x, y, z};
            ctrlnet << x << " " << y << " " << z << "\n";
        }
    }

    BSplineSurface surf(ub, vb, ctrl);

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

    std::cout << "Wrote plots/data/surface.csv and tsurface.xyz\n";
    std::cout << "Demo complete.  Wrote " << ptCount << " points.\n";

    auto Su = surf.derivative_u(0.5, 0.5);
    auto Sv = surf.derivative_v(0.5, 0.5);
    auto N = surf.normal(0.5, 0.5);

    std::cout << std::endl
              << "Derivatives at (u,v) = (0.5,0.5):" << std::endl;
    std::cout << "Su = " << Su[0] << "," << Su[1] << "," << Su[2] << "\n";
    std::cout << "Sv = " << Sv[0] << "," << Sv[1] << "," << Sv[2] << "\n";
    std::cout << "N  = " << N[0] << "," << N[1] << "," << N[2] << "\n";

    auto c = surf.curvature(0.5, 0.5);
    std::cout << std::endl
              << "Curvatures at (u,v) = (0.5,0.5):" << std::endl;
    std::cout << "K = " << c.K << " H = " << c.H << " k1=" << c.k1 << " k2=" << c.k2 << "\n";

    std::cout << std::endl
              << "Curvatures at (u,v) = (0.5,0.5):" << std::endl;
    std::cout << "K = " << c.K << " H = " << c.H << " k1=" << c.k1 << " k2=" << c.k2 << "\n";

    std::cout << std::endl
              << "Surface info dump:" << std::endl;
    surf.DumpInfo();

    std::cout << std::endl
              << "=== Running Closest Point Test ===" << std::endl;
    test_closest_point(surf);

    return 0;
}
