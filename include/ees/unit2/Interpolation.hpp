#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include "ees/core/NumericalResult.hpp"
#include <vector>
#include <string>

namespace ees::unit2 {

struct DifferenceTableResult {
    core::Vector x;
    core::Vector y;
    core::Matrix diffTable; // Column 0 = y, Column 1 = Delta y, Column 2 = Delta^2 y, etc.
    bool isEquispaced{false};
    double h{0.0};
    std::string toString() const;
};

struct InterpolationResult {
    double evaluatedX{0.0};
    double interpolatedY{0.0};
    double estimatedError{0.0};
    std::string method;
    std::string polynomialString; // e.g. "P(x) = 2.0 + 3.0*(x - 1) + ..."
    std::vector<double> coefficients;
    std::string summary() const;
};

class InterpolationLab {
public:
    // Finite Differences
    static DifferenceTableResult computeForwardDifferenceTable(
        const core::Vector& x,
        const core::Vector& y
    );

    static DifferenceTableResult computeBackwardDifferenceTable(
        const core::Vector& x,
        const core::Vector& y
    );

    static DifferenceTableResult computeDividedDifferenceTable(
        const core::Vector& x,
        const core::Vector& y
    );

    // Gregory-Newton Forward Interpolation
    // Formula: y(x) = y0 + u*Delta y0 + u(u-1)/2! * Delta^2 y0 + ... where u = (x - x0)/h
    static InterpolationResult gregoryNewtonForward(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );

    // Gregory-Newton Backward Interpolation
    // Formula: y(x) = yn + v*Nabla yn + v(v+1)/2! * Nabla^2 yn + ... where v = (x - xn)/h
    static InterpolationResult gregoryNewtonBackward(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );

    // Lagrange Interpolation
    // Formula: P(x) = sum_{i=0}^n y_i * L_i(x), where L_i(x) = prod_{j!=i} (x - x_j)/(x_i - x_j)
    static InterpolationResult lagrange(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );

    // Newton Divided Difference Interpolation
    // Formula: P(x) = f[x0] + f[x0,x1]*(x-x0) + f[x0,x1,x2]*(x-x0)*(x-x1) + ...
    static InterpolationResult newtonDividedDifference(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );

    // Compare all interpolation methods on dataset
    static std::string compareInterpolation(
        const core::Vector& x,
        const core::Vector& y,
        double xVal
    );
};

} // namespace ees::unit2

