#include "test_closest_point.hpp"
#include <random>
#include <iostream>
#include <iomanip>
#include "bspline_surface.hpp"

struct Stats
{
    int total = 0;
    int failures = 0;
    double maxDistError = 0.0;
    double sumDistError = 0.0;
    double maxParamError = 0.0;
    double sumParamError = 0.0;
};

void test_closest_point(BSplineSurface &surf)
{
    Stats stats;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.05); // 5cm noise, adjust as needed

    const int N = 2000; // number of random tests

    for (int k = 0; k < N; ++k)
    {
        stats.total++;

        // 1. Pick a random true parameter
        double u_true = uni(rng);
        double v_true = uni(rng);

        // 2. Evaluate true surface point
        auto P = surf.evaluate(u_true, v_true);

        // 3. Perturb it to create a test point Q
        std::array<double, 3> Q = {
            P[0] + noise(rng),
            P[1] + noise(rng),
            P[2] + noise(rng)};

        // 4. Initial guess (random or center)
        double u0 = uni(rng);
        double v0 = uni(rng);

        // 5. Run LM closest-point
        auto uv = surf.closest_point_LM(Q, u0, v0, 30, 1e-10);
        double u_est = uv.first;
        double v_est = uv.second;

        // 6. Evaluate recovered point
        auto Pest = surf.evaluate(u_est, v_est);

        // 7. Distance from Q to recovered closest point
        double dxQ = Pest[0] - Q[0];
        double dyQ = Pest[1] - Q[1];
        double dzQ = Pest[2] - Q[2];
        double distToQ = std::sqrt(dxQ * dxQ + dyQ * dyQ + dzQ * dzQ);

        // 8. For curiosity: distance from original point P to Q
        double dxP = P[0] - Q[0];
        double dyP = P[1] - Q[1];
        double dzP = P[2] - Q[2];
        double distPtoQ = std::sqrt(dxP * dxP + dyP * dyP + dzP * dzP);

        // 9. Parameter error (still fine as a diagnostic)
        double du = std::abs(u_est - u_true);
        double dv = std::abs(v_est - v_true);
        double paramErr = std::max(du, dv);

        // 10. Update stats using distToQ
        stats.sumDistError += distToQ;
        stats.maxDistError = std::max(stats.maxDistError, distToQ);
        stats.sumParamError += paramErr;
        stats.maxParamError = std::max(stats.maxParamError, paramErr);

        // 11. Failure detection: did we *reduce* the distance?
        if (distToQ > distPtoQ + 1e-6) // projection made it worse
            stats.failures++;
    }

    // Print summary
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Total tests: " << stats.total << "\n";
    std::cout << "Failures:    " << stats.failures << "\n";
    std::cout << "Mean dist error: " << stats.sumDistError / stats.total << "\n";
    std::cout << "Max dist error:  " << stats.maxDistError << "\n";
    std::cout << "Mean param error: " << stats.sumParamError / stats.total << "\n";
    std::cout << "Max param error:  " << stats.maxParamError << "\n";
}