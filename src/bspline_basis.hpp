#pragma once
#include <vector>

class BSplineBasis
{
public:
    BSplineBasis(int degree, std::vector<double> knots);

    // Evaluate N_i^p(u)
    double evaluate(int i, int p, double u) const;

    // Derivatives
    double derivative(int i, int p, double u) const;
    double derivative(int i, double u) const { return derivative(i, degree_, u); }
    double second_derivative(int i, int p, double u) const;
    double second_derivative(int i, double u) const { return second_derivative(i, degree_, u); }

    int degree() const { return degree_; }
    const std::vector<double> &knots() const { return knots_; }

    // Convenience: evaluate using stored degree
    double evaluate(int i, double u) const { return evaluate(i, degree_, u); }

    // Number of basis functions == number of control points this basis expects
    int numBasis() const { return static_cast<int>(knots_.size()) - degree_ - 1; }

    // Get extents of basis functions to [umin, umax]
    std::pair<double, double> getExtents() const;

    static void check_basis_consistency(const BSplineBasis &B, int numCtrl);

    void DumpInfo(const char *msg) const;
    void GeneratePlot(const char *filename) const;

    static constexpr double DOMAIN_PADDING_FACTOR = 1e-6;

private:
    int degree_;
    std::vector<double> knots_;
};
