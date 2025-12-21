#include "bspline_basis.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

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
