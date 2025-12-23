#include <random>
#include <iostream>
#include <iomanip>
#include "bspline_surface_builder.hpp"
#include <chrono>

struct Stats
{
    int total = 0;
    int failures = 0;
    int edgeCases = 0;
    int stagnations = 0;
    int divergences = 0;
    int maxiterations = 0;
    int iterations = 0;
    int numValid = 0;
    double maxDistError = 0.0;
    double sumDistError = 0.0;
    double maxParamError = 0.0;
    double sumParamError = 0.0;
    double maxObjective = 0.0;
    double sumObjective = 0.0;
};

void test_closest_point(BSplineSurface &surf)
{
    Stats stats;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.05); // 5cm noise, adjust as needed

    const int N = 100000; // number of random tests
    auto start = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        stats.total++;

        if (stats.total % 1000 == 0)
        {
            auto end = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration<double>(end - start).count();
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "Test " << stats.total << " of " << N << " (" << elapsed << "seconds)" << std::endl;
        }

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
        auto result = surf.closest_point_global(Q, 4, 0.5, 0.5, 100, 1e-8, 100, 1e-8);
        double u_est = result.u;
        double v_est = result.v;
        stats.iterations += result.iterations;
        if (result.status == ClosestPointResult::Status::Boundary)
        {
            stats.edgeCases++;
        }
        else if (result.status != ClosestPointResult::Status::Success)
        {
            if (result.status == ClosestPointResult::Status::MaxIterations)
            {
                stats.maxiterations++;
            }
            else if (result.status == ClosestPointResult::Status::Stagnation)
            {
                stats.stagnations++;
                stats.maxObjective = std::max(stats.maxObjective, result.bestObjective);
                stats.sumObjective += result.bestObjective;
                stats.numValid++;
            }
            else if (result.status == ClosestPointResult::Status::Divergence)
            {
                stats.divergences++;
            }
            else
            {
                stats.maxObjective = std::max(stats.maxObjective, result.bestObjective);
                stats.sumObjective += result.bestObjective;
                stats.numValid++;
            }
        }

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
    std::cout << "Total tests     : " << stats.total << "\n";
    std::cout << "Failures        : " << stats.failures << "\n";
    std::cout << "Edge cases      : " << stats.edgeCases << "\n";
    std::cout << "Stagnations     : " << stats.stagnations << "\n";
    std::cout << "Divergences     : " << stats.divergences << "\n";
    std::cout << "Max iterations  : " << stats.maxiterations << "\n";
    std::cout << "Avg Iterations  : " << double(stats.iterations) / stats.total << "\n";
    std::cout << "Mean dist error : " << stats.sumDistError / stats.total << "\n";
    std::cout << "Max dist error  : " << stats.maxDistError << "\n";
    std::cout << "Mean param error: " << stats.sumParamError / stats.total << "\n";
    std::cout << "Max param error : " << stats.maxParamError << "\n";

    if (stats.numValid > 0)
    {
        std::cout << "Mean objective  : " << stats.sumObjective / stats.numValid << "\n";
        std::cout << "Max objective   : " << stats.maxObjective << "\n";
    }
}

int main()
{
    std::cout << "=== Testing closest point ===\n";
    auto surf = SurfaceBuilder::Bowl();
    test_closest_point(surf);
    return (0);
}
