#include "bspline_basis.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>

BSplineBasis::BSplineBasis(
    int degree,
    std::vector<double> knots,
    bool isPeriodic /* = false */)
    : degree_(degree), knots_(std::move(knots)), isPeriodic_(isPeriodic)
{
    nBasis_ = static_cast<int>(knots_.size()) - degree_ - 1;

    if (nBasis_ <= 0)
        throw std::runtime_error("BSplineBasis: invalid basis");
}

double BSplineBasis::evaluate(int i, int p, double u) const
{
    // Number of basis functions for this basis (fixed for all p)
    const int maxI = nBasis_ - 1;

    // Basic safety: i must be valid for N_i^p
    if (isPeriodic_)
    {
        i = wrapIndex(i);
    }
    else
    {
        if (i < 0 || i > maxI)
            return 0.0;
    }

    const auto &U = knots_;

    if (p == 0)
    {
        // Standard interior case: [U[i], U[i+1))
        if (u >= U[i] && u < U[i + 1])
            return 1.0;

        // Special case: right endpoint u == last knot -> assign to last basis
        if (!isPeriodic_)
        {
            if (u == U.back() && i == maxI)
                return 1.0;
        }

        return 0.0;
    }

    double left = 0.0;
    double right = 0.0;

    const double denom1 = U[i + p] - U[i];
    const double denom2 = U[i + p + 1] - U[i + 1];

    if (denom1 > EPSILON)
        left = (u - U[i]) / denom1 * evaluate(i, p - 1, u);

    if (denom2 > EPSILON)
        right = (U[i + p + 1] - u) / denom2 * evaluate((isPeriodic_ ? wrapIndex(i + 1) : i + 1), p - 1, u);

    return left + right;
}

double BSplineBasis::derivative(int i, int p, double u) const
{
    // Number of basis functions for this basis (fixed for all p)
    const int maxI = nBasis_ - 1;

    // Same bounds check as evaluate()
    if (isPeriodic_)
    {
        i = wrapIndex(i);
    }
    else
    {
        if (i < 0 || i > maxI)
            return 0.0;
    }

    const auto &U = knots_;

    if (p == 0)
        return 0.0; // degree 0 basis has zero derivative

    double term1 = 0.0;
    double term2 = 0.0;

    double denom1 = U[i + p] - U[i];
    double denom2 = U[i + p + 1] - U[i + 1];

    if (denom1 > EPSILON)
        term1 = (p / denom1) * evaluate(i, p - 1, u);

    if (denom2 > EPSILON)
        term2 = (p / denom2) * evaluate((isPeriodic_ ? wrapIndex(i + 1) : i + 1), p - 1, u);

    return term1 - term2;
}

double BSplineBasis::second_derivative(int i, int p, double u) const
{
    // Number of basis functions for this basis (fixed for all p)
    const int maxI = nBasis_ - 1;

    if (isPeriodic_)
    {
        i = wrapIndex(i);
    }
    else
    {
        if (i < 0 || i > maxI)
            return 0.0;
    }

    const auto &U = knots_;

    if (p == 0)
        return 0.0; // piecewise constant -> second derivative zero

    double denom1 = U[i + p] - U[i];
    double denom2 = U[i + p + 1] - U[i + 1];

    double term1 = 0.0;
    double term2 = 0.0;

    if (denom1 > EPSILON)
        term1 = (p / denom1) * derivative(i, p - 1, u);

    if (denom2 > EPSILON)
        term2 = (p / denom2) * derivative((isPeriodic_ ? wrapIndex(i + 1) : i + 1), p - 1, u);

    return term1 - term2;
}

std::pair<double, double> BSplineBasis::getExtents() const
{
    if (degree_ < 1)
        throw std::runtime_error("BSplineBasis: Degree must be >= 1 for geometric curves");

    if (knots_.size() < 2 * degree_ + 2)
        throw std::runtime_error("BSplineBasis: Invalid knot vector");

    int iMin = degree_;
    int iMax = nBasis_;

    if (iMax <= iMin)
        throw std::runtime_error("BSplineBasis: degenerate knot span");

    double tMin = knots_[iMin];
    double tMax = knots_[iMax];

    if (tMin < tMax)
        return {tMin, tMax};

    throw std::runtime_error("BSplineBasis: invalid domain");
}

void BSplineBasis::DumpInfo(const char *msg /* = nullptr */) const
{
    const int nKnots = static_cast<int>(knots_.size());

    std::cout << ((msg != nullptr) ? msg : "") << "B-Spline Basis Info:\n";
    std::cout << " Degree: " << degree_ << "\n";
    std::cout << " Num Basis Functions: " << nBasis_ << "\n";
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

    auto u_extents = getExtents();
    double u_delta = (u_extents.second - u_extents.first) / 100;

    for (int i = 0; i <= 100; i++)
    {
        double u = u_extents.first + i * u_delta;
        u = std::clamp(u, u_extents.first, u_extents.second);

        outFile << u;

        for (int i = 0; i < numBasis(); ++i)
            outFile << "," << evaluate(i, u);

        outFile << std::endl;
    }

    outFile.close();

    std::cout << "Wrote basis plot file: " << filename << std::endl;
}
