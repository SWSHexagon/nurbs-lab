#include <random>
#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>
#include <numbers>
#include <fstream>
#include <exception>
#include "bspline_curve_builder.hpp"
#include "line_curve.hpp"
#include "circle_curve.hpp"
#include "nurbs_curve_builder.hpp"

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

    double sumOffset = 0.0;

    double sumStagnationGradients = 0.0;
};

void test_closest_point_curve(ParametricCurve &curve)
{
    CurveStats stats;

    auto domain = curve.domain();
    double t_start = domain.first;
    double t_end = domain.second;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> uni(t_start, t_end);
    std::normal_distribution<double> noise(0.0, 0.02); // 2cm noise
    std::uniform_real_distribution<double> radialOffset(0.8, 1.2);

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
        auto Q = scale(P, radialOffset(rng));

        // 4. Random initial guess
        // double t0 = t_true + 0.2 * uni(rng);
        double t0 = uni(rng);

        // 5. Run global closest-point
        auto result = curve.closest_point_LM(Q, t0, 100, 1e-8);

        stats.iterations += result.iterations;

        if (result.status == ClosestPointResult::Status::MaxIterations)
            stats.maxiterations++;
        else if (result.status == ClosestPointResult::Status::Stagnation)
        {
            stats.stagnations++;
            stats.maxObjective = std::max(stats.maxObjective, result.bestObjective);
            stats.sumObjective += result.bestObjective;
            stats.numValid++;
            stats.sumStagnationGradients += result.gradNorm;
        }
        else if (result.status == ClosestPointResult::Status::Divergence)
            stats.divergences++;

        // 6. Evaluate recovered point
        auto Pest = result.point3D;

        // 7. Distance from resolved point to Q
        double dx = Pest[0] - P[0];
        double dy = Pest[1] - P[1];
        double dz = Pest[2] - P[2];
        double distErr = std::sqrt(dx * dx + dy * dy + dz * dz);

        dx = Pest[0] - Q[0];
        dy = Pest[1] - Q[1];
        dz = Pest[2] - Q[2];
        double ptOffset = std::sqrt(dx * dx + dy * dy + dz * dz);

        stats.sumOffset += ptOffset;

        // 9. Parameter error
        double raw = std::abs(result.u - t_true);
        double period = domain.second - domain.first;
        double paramErr = std::min(raw, period - raw);

        // 10. Update stats
        stats.sumDistError += distErr;
        stats.maxDistError = std::max(stats.maxDistError, distErr);

        stats.sumParamError += paramErr;
        stats.maxParamError = std::max(stats.maxParamError, paramErr);

        // 11. Failure detection
        if (distErr > 1e-3)
            stats.failures++;
    }

    // Summary
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Range             : (" << t_start << ", " << t_end << ")" << std::endl;
    std::cout << "Tolerance         : " << 1e-8 << std::endl;
    std::cout << "Total tests       : " << stats.total << "\n";
    std::cout << "Failures          : " << stats.failures << "\n";
    std::cout << "Stagnations       : " << stats.stagnations << "\n";
    std::cout << "Divergences       : " << stats.divergences << "\n";
    std::cout << "Max iterations    : " << stats.maxiterations << "\n";
    std::cout << "Avg Iterations    : " << double(stats.iterations) / stats.total << "\n";
    std::cout << "Mean point error  : " << stats.sumDistError / stats.total << "\n";
    std::cout << "Max point error   : " << stats.maxDistError << "\n";
    std::cout << "Mean point offset : " << stats.sumOffset / stats.total << "\n";
    std::cout << "Mean param error  : " << stats.sumParamError / stats.total << "\n";
    std::cout << "Max param error   : " << stats.maxParamError << "\n";

    if (stats.numValid > 0)
    {
        std::cout << "Mean objective  : " << stats.sumObjective / stats.numValid << "\n";
        std::cout << "Max objective   : " << stats.maxObjective << "\n";
    }

    if (stats.stagnations > 0)
    {
        std::cout << "Avg stag grads  : " << stats.sumStagnationGradients / stats.stagnations << "\n";
    }

    auto end = std::chrono::steady_clock::now();
    double totalSeconds = std::chrono::duration<double>(end - start).count();
    std::cout << "Total elapsed   : " << totalSeconds << " seconds\n";
    std::cout << "Avg per solve   : " << (totalSeconds / N) * 1e6 << " microseconds\n";
    std::cout << "Solves per sec  : " << N / totalSeconds << "\n";
}

void generate_data_file(ParametricCurve &curve, const char *fname)
{
    std::ofstream outFile(fname);

    if (!outFile)
    {
        std::cout << "Unable to open file: " << fname << std::endl;
        return;
    }

    const int N = 1000;

    auto domain = curve.domain();
    double t_min = domain.first;
    double t_max = domain.second;
    double t_delta = (t_max - t_min) / N;

    for (int i = 0; i < N; i++)
    {
        double t = t_min + i * t_delta;
        auto p = curve.evaluate(t);
        outFile << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }

    outFile.close();
}

int main()
{
    std::cout << "=== Testing curve closest point ===\n";

    try
    {
        // auto curve = CurveBuilder::Stress(); // or any analytic test curve
        // LineCurve curve({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
        CircleCurve curve({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0);
        // NURBSCurve curve = NURBSCurveBuilder::MakeNURBSCircle({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0);
        // NURBSCurve curve = NURBSCurveBuilder::MakeNURBSPeriodicCircle(1.0);
        // NURBSCurve curve = NURBSCurveBuilder::MakeNURBSCubicCircle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.1}, 1.0);

        test_closest_point_curve(curve);
        generate_data_file(curve, "C:\\labs\\nurbs-lab\\plots\\data\\curve.xyz");

        auto domain = curve.domain();
        double domainDelta = (domain.second - domain.first) / 8;

        std::cout << "Derivatives at 45 degree intervals: " << std::endl;
        for (int i = 0; i < 9; i++)
        {
            double t = domain.first + double(i) * domainDelta;

            auto d1 = curve.derivative(t);
            auto d2 = curve.second_derivative(t);

            std::cout << "t: " << t << ", " << "d1: (" << d1[0] << ", " << d1[1] << ", " << d1[2] << "), d2: (" << d2[0] << ", " << d2[1] << ", " << d2[2] << ")" << std::endl;
        }

        std::cout << "Domain          : [" << domain.first << ", " << domain.second << "]" << std::endl;

        curve.DumpInfo();
    }
    catch (const std::exception &e)
    {
        std::cout << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
