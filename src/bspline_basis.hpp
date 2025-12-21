#pragma once
#include <vector>

class BSplineBasis
{
public:
    BSplineBasis(int degree, std::vector<double> knots);

    // Evaluate N_i^p(u)
    double evaluate(int i, int p, double u) const;

    int degree() const { return degree_; }
    const std::vector<double> &knots() const { return knots_; }

    // Convenience: evaluate using stored degree
    double evaluate(int i, double u) const { return evaluate(i, degree_, u); }

    // Number of basis functions == number of control points this basis expects
    int numBasis() const { return static_cast<int>(knots_.size()) - degree_ - 1; }

private:
    int degree_;
    std::vector<double> knots_;
};