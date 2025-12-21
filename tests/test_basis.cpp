#include "bspline_basis.hpp"
#include <cassert>
#include <iostream>

int main()
{
    BSplineBasis basis(2, {0, 0, 0, 1, 2, 3, 3, 3});
    double val = basis.evaluate(2, 2, 1.0);
    assert(val >= 0.0 && val <= 1.0);
    std::cout << "Basis test passed\n";

    return (0);
}