#pragma once

#include <array>

struct ClosestPointResult
{
    double u = 0.0;
    double v = 0.0;

    std::array<double, 3> point3D; // final closest point on surface

    double distance = 0.0; // final |S(u,v) - P|
    double gradNorm = 0.0; // ||∇(0.5||F||²)||

    int iterations = 0; // final number of iterations

    enum class Status
    {
        Success,      // Interior convergence
        Boundary,     // Converged on boundary
        Stagnation,   // Step size too small with large gradient norm
        Divergence,   // lambda too large or Hessian unusable
        MaxIterations // Loop exhaused without convergence
    } status = Status::MaxIterations;

    bool onBoundary = false; // Whether the final point is on the boundary
};
