#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include <vector>
#include <string>

namespace ees::unit3 {

struct CubicSplineSegment {
    double x_i{0.0};
    double a{0.0}; // constant
    double b{0.0}; // linear
    double c{0.0}; // quadratic
    double d{0.0}; // cubic: S_i(x) = a + b*(x-xi) + c*(x-xi)^2 + d*(x-xi)^3
};

struct SplineResult {
    double evaluatedX{0.0};
    double interpolatedY{0.0};
    std::string splineType;
    std::vector<CubicSplineSegment> segments; // Populated for cubic splines
    std::string summary() const;
};

class SplineLab {
public:
    // First-degree Spline (Linear Spline)
    // S_i(x) = y_i + (y_{i+1} - y_i)/(x_{i+1} - x_i) * (x - x_i)
    static SplineResult linearSpline(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );

    // Second-degree Spline (Quadratic Spline)
    // Continuous 0th and 1st derivatives across interior knots
    static SplineResult quadraticSpline(
        const core::Vector& x,
        const core::Vector& y,
        double xVal,
        double initialM0 = 0.0 // Boundary condition S'_0(x_0) = initialM0
    );

    // Natural Cubic Spline
    // S''(x_0) = 0 and S''(x_n) = 0 (Natural / Free boundary conditions)
    // Solves tridiagonal linear system for second derivatives using Thomas algorithm
    static SplineResult naturalCubicSpline(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );

    // Natural Cubic Spline evaluation across multiple points (for plotting curves)
    static core::Vector evaluateNaturalCubicSpline(
        const core::Vector& xKnots,
        const core::Vector& yKnots,
        const core::Vector& evalPoints
    );

    // Uniform B-Spline curve generation (Cox-de Boor algorithm)
    // degree p (typically cubic p=3)
    static core::Vector evaluateBSpline(
        const core::Vector& controlPointsX,
        const core::Vector& controlPointsY,
        int numSamples = 100,
        int degree = 3
    );

    // Compare all spline representations against original data
    static std::string compareSplines(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );
};

} // namespace ees::unit3

