#pragma once
#include <vector>
#include <limits>

class BSplineBasis
{
public:
    BSplineBasis(
        int degree,
        std::vector<double> knots,
        bool isPeriodic = false);

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
    inline int numBasis() const { return nBasis_; }

    // Get extents of basis functions to [umin, umax]
    std::pair<double, double> getExtents() const;

    void DumpInfo(const char *msg = nullptr) const;
    void GeneratePlot(const char *filename) const;

    static constexpr double EPSILON = 100 * std::numeric_limits<double>::epsilon();

private:
    int nBasis_ = 0;
    int degree_ = 2;
    bool isPeriodic_ = false;
    std::vector<double> knots_;

    inline int wrapIndex(int i) const
    {
        return (i % nBasis_ + nBasis_) % nBasis_;
    }
};
