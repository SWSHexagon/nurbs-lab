#pragma once

#include "parametric_surface.hpp"
#include "bspline_basis.hpp"
#include "surface_differential.hpp"
#include "surface_curvature.hpp"
#include "surface_frame.hpp"
#include "bspline_curve.hpp"
#include <vector>
#include <array>
#include <utility>
#include <tuple>

class BSplineSurface : public ParametricSurface
{
public:
    BSplineSurface(
        BSplineBasis u_basis,
        BSplineBasis v_basis,
        std::vector<std::vector<std::array<double, 3>>> control_points);

    virtual std::array<double, 3> evaluate(double u, double v) const;

    virtual std::pair<std::array<double, 3>, std::array<double, 3>> derivatives(
        double u, double v) const { return {derivative_u(u, v), derivative_v(u, v)}; }

    virtual std::array<double, 4> domain() const;

    virtual std::array<double, 3> normal(double u, double v, bool unitize = true) const;

    SurfaceDifferential differential(double u, double v) const;

    SurfaceCurvature curvature(double u, double v) const;

    SurfaceFrame frame(double u, double v) const;

    virtual const BSplineCurve &boundary_u_min() const { return boundary_u_min_; }
    virtual const BSplineCurve &boundary_u_max() const { return boundary_u_max_; }
    virtual const BSplineCurve &boundary_v_min() const { return boundary_v_min_; }
    virtual const BSplineCurve &boundary_v_max() const { return boundary_v_max_; }

    void DumpInfo() const;

private:
    BSplineBasis u_basis_;
    BSplineBasis v_basis_;
    std::vector<std::vector<std::array<double, 3>>> ctrl_;

    // Boundary curves
    BSplineCurve boundary_u_min_{BSplineBasis(0, {0, 0}), {}};
    BSplineCurve boundary_u_max_{BSplineBasis(0, {0, 0}), {}};
    BSplineCurve boundary_v_min_{BSplineBasis(0, {0, 0}), {}};
    BSplineCurve boundary_v_max_{BSplineBasis(0, {0, 0}), {}};

    // Helper function to initialize boundary curves from control points
    void build_boundary_curves();

    // First derivatives of the surface
    std::array<double, 3> derivative_u(double u, double v) const;
    std::array<double, 3> derivative_v(double u, double v) const;

    // Second derivatives of the surface
    std::array<double, 3> second_derivative_uu(double u, double v) const;
    std::array<double, 3> second_derivative_vv(double u, double v) const;
    std::array<double, 3> second_derivative_uv(double u, double v) const;

    static std::array<double, 2> principal_direction_uv(
        const SurfaceDifferential &d,
        double k);

    static std::array<double, 3> combine_uv_to_3d(
        const SurfaceDifferential &d,
        const std::array<double, 2> &uv);
};