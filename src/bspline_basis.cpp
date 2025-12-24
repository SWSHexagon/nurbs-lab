#include "bspline_basis.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>

BSplineBasis::BSplineBasis(int degree, std::vector<double> knots)
    : degree_(degree), knots_(std::move(knots))
{
    assert(degree_ >= 0);
    assert(knots_.size() >= static_cast<size_t>(degree_ + 2));
}

double BSplineBasis::evaluate(int i, int p, double u) const
{
    const int nKnots = static_cast<int>(knots_.size());
    const auto &U = knots_;

    // Number of basis functions for this basis (fixed for all p)
    const int nBasis = nKnots - degree_ - 1;
    const int maxI = nBasis - 1;

    // Basic safety: i must be valid for N_i^p
    if (i < 0 || i > maxI)
        return 0.0;

    if (p == 0)
    {
        // Standard interior case: [U[i], U[i+1))
        if (u >= U[i] && u < U[i + 1])
            return 1.0;

        // Special case: right endpoint u == last knot -> assign to last basis
        if (u == U.back() && i == maxI)
            return 1.0;

        return 0.0;
    }

    double left = 0.0;
    double right = 0.0;

    const double denom1 = U[i + p] - U[i];
    const double denom2 = U[i + p + 1] - U[i + 1];

    if (denom1 > 1e-14)
        left = (u - U[i]) / denom1 * evaluate(i, p - 1, u);

    if (denom2 > 1e-14)
        right = (U[i + p + 1] - u) / denom2 * evaluate(i + 1, p - 1, u);

    return left + right;
}

double BSplineBasis::derivative(int i, int p, double u) const
{
    const int nKnots = knots_.size();
    const int nBasis = nKnots - degree_ - 1;
    const int maxI = nBasis - 1;

    // Same bounds check as evaluate()
    if (i < 0 || i > maxI)
        return 0.0;

    const auto &U = knots_;

    if (p == 0)
        return 0.0; // degree 0 basis has zero derivative

    double term1 = 0.0;
    double term2 = 0.0;

    double denom1 = U[i + p] - U[i];
    double denom2 = U[i + p + 1] - U[i + 1];

    if (denom1 > 1e-14)
        term1 = (p / denom1) * evaluate(i, p - 1, u);

    if (denom2 > 1e-14)
        term2 = (p / denom2) * evaluate(i + 1, p - 1, u);

    return term1 - term2;
}

double BSplineBasis::second_derivative(int i, int p, double u) const
{
    const int nKnots = static_cast<int>(knots_.size());
    const int nBasis = nKnots - degree_ - 1;
    const int maxI = nBasis - 1;

    if (i < 0 || i > maxI)
        return 0.0;

    const auto &U = knots_;

    if (p == 0)
        return 0.0; // piecewise constant -> second derivative zero

    double denom1 = U[i + p] - U[i];
    double denom2 = U[i + p + 1] - U[i + 1];

    double term1 = 0.0;
    double term2 = 0.0;

    if (denom1 > 1e-14)
        term1 = (p / denom1) * derivative(i, p - 1, u);

    if (denom2 > 1e-14)
        term2 = (p / denom2) * derivative(i + 1, p - 1, u);

    return term1 - term2;
}

std::pair<double, double> BSplineBasis::getExtents() const
{
    if (knots_.size() < 2 * degree_ + 2)
        throw std::runtime_error("Invalid knot vector");

    int iMin = degree_;
    int iMax = static_cast<int>(knots_.size()) - degree_ - 1;

    if ((iMin >= 0) && (iMax >= iMin))
    {
        return std::pair<double, double>(knots_[iMin], knots_[iMax]);
    }

    return std::pair<double, double>();
}

void BSplineBasis::check_basis_consistency(const BSplineBasis &B, int numCtrl)
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

void BSplineBasis::DumpInfo(const char *msg) const
{
    const int nBasis = numBasis();
    const int nKnots = static_cast<int>(knots_.size());

    std::cout << msg << "B-Spline Basis Info:\n";
    std::cout << " Degree: " << degree_ << "\n";
    std::cout << " Num Basis Functions: " << nBasis << "\n";
    std::cout << " Knot Vector (" << nKnots << " knots): ";
    for (double k : knots_)
        std::cout << k << " ";
    std::cout << "\n";
}

void BSplineBasis::GeneratePlot(const char *filename) const
{
    std::ofstream outFile(filename);

    if (!outFile)
    {
        std::cerr << "Failed to open basis plot file:" << filename << std::endl;
        return;
    }

    for (int i = 0; i <= 100; i++)
    {
        double u = 0.01 * double(i);

        if (u > 1.0)
            u = 1.0;

        outFile << u;

        for (int i = 0; i < numBasis(); ++i)
            outFile << "," << evaluate(i, u);

        outFile << std::endl;
    }

    outFile.close();

    std::cout << "Wrote basis plot file: " << filename << std::endl;
}
