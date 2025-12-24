#include <random>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "bspline_curve_builder.hpp"
#include "line_curve.hpp"
#include "circle_curve.hpp"

struct CurveStats
{
    int total = 0;
    int failures = 0;
    int stagnations = 0;
    int divergences = 0;
    int maxiterations = 0;
    int numValid = 0;

    int iterations = 0;

    double maxDistError = 0.0;
    double sumDistError = 0.0;

    double maxParamError = 0.0;
    double sumParamError = 0.0;

    double maxObjective = 0.0;
    double sumObjective = 0.0;
};

void test_closest_point_curve(ParametricCurve &curve)
{
    CurveStats stats;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.02); // 2cm noise

    const int N = 100000;
    auto start = std::chrono::steady_clock::now();

    for (int k = 0; k < N; ++k)
    {
        stats.total++;

        if (stats.total % 1000 == 0)
        {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start).count();
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "Test " << stats.total << " of " << N
                      << " (" << elapsed << " seconds)\n";
        }

        // 1. Random true parameter
        double t_true = uni(rng);

        // 2. True point on curve
        auto P = curve.evaluate(t_true);

        // 3. Perturb to create noisy query Q
        std::array<double, 3> Q = {
            P[0] + noise(rng),
            P[1] + noise(rng),
            P[2] + noise(rng)};

        // 4. Random initial guess
        double t0 = uni(rng);

        // 5. Run global closest-point
        auto result = curve.closest_point_LM(Q, 0.5, 100, 1e-8);

        stats.iterations += result.iterations;

        if (result.status == ClosestPointResult::Status::MaxIterations)
            stats.maxiterations++;
        else if (result.status == ClosestPointResult::Status::Stagnation)
        {
            stats.stagnations++;
            stats.maxObjective = std::max(stats.maxObjective, result.bestObjective);
            stats.sumObjective += result.bestObjective;
            stats.numValid++;
        }
        else if (result.status == ClosestPointResult::Status::Divergence)
            stats.divergences++;

        // 6. Evaluate recovered point
        auto Pest = curve.evaluate(result.u);

        // 7. Distance from Q to recovered closest point
        double dx = Pest[0] - Q[0];
        double dy = Pest[1] - Q[1];
        double dz = Pest[2] - Q[2];
        double distToQ = std::sqrt(dx * dx + dy * dy + dz * dz);

        // 8. Distance from true point to Q
        double dxP = P[0] - Q[0];
        double dyP = P[1] - Q[1];
        double dzP = P[2] - Q[2];
        double distPtoQ = std::sqrt(dxP * dxP + dyP * dyP + dzP * dzP);

        // 9. Parameter error
        double paramErr = std::abs(result.u - t_true);

        // 10. Update stats
        stats.sumDistError += distToQ;
        stats.maxDistError = std::max(stats.maxDistError, distToQ);

        stats.sumParamError += paramErr;
        stats.maxParamError = std::max(stats.maxParamError, paramErr);

        // 11. Failure detection
        if (distToQ > distPtoQ + 1e-6)
            stats.failures++;
    }

    // Summary
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Total tests     : " << stats.total << "\n";
    std::cout << "Failures        : " << stats.failures << "\n";
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

    auto end = std::chrono::steady_clock::now();
    double totalSeconds = std::chrono::duration<double>(end - start).count();
    std::cout << "Total elapsed   : " << totalSeconds << " seconds\n";
    std::cout << "Avg per solve   : " << (totalSeconds / N) * 1e6 << " microseconds\n";
    std::cout << "Solves per sec  : " << N / totalSeconds << "\n";
}

int main()
{
    std::cout << "=== Testing curve closest point ===\n";
    // auto curve = CurveBuilder::Stress(); // or any analytic test curve
    // LineCurve curve({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    CircleCurve curve({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0);
    test_closest_point_curve(curve);
    return 0;
}
