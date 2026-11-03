#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include "ees/core/NumericalResult.hpp"
#include "ees/core/ExpressionParser.hpp"
#include <vector>
#include <string>

namespace ees::unit2 {

struct IntegrationResult {
    double value{0.0};
    double trueValue{0.0};
    double absError{0.0};
    double relError{0.0};
    int subintervals{1};
    int functionEvaluations{0};
    std::string method;
    core::Matrix rombergTable; // Only populated for Romberg
    std::string summary() const;
};

class IntegrationLab {
public:
    // Trapezoidal Rule: composite over n subintervals
    // Formula: (h/2) * [f(a) + 2*sum(f(x_i)) + f(b)], Error: O(h^2)
    static IntegrationResult trapezoidal(
        core::RealFunction f,
        double a,
        double b,
        int subintervals = 100,
        double trueVal = 0.0
    );

    // Simpson's 1/3 Rule: composite over n subintervals (n must be even)
    // Formula: (h/3) * [f(a) + 4*sum_odd + 2*sum_even + f(b)], Error: O(h^4)
    static IntegrationResult simpson13(
        core::RealFunction f,
        double a,
        double b,
        int subintervals = 100,
        double trueVal = 0.0
    );

    // Simpson's 3/8 Rule: composite over n subintervals (n must be multiple of 3)
    // Formula: (3h/8) * [f(a) + 3*f(x1) + 3*f(x2) + 2*f(x3) + ... + f(b)], Error: O(h^4)
    static IntegrationResult simpson38(
        core::RealFunction f,
        double a,
        double b,
        int subintervals = 99,
        double trueVal = 0.0
    );

    // Romberg Algorithm: Richardson extrapolation applied to Trapezoidal rule
    // R(k, j) = R(k, j-1) + [R(k, j-1) - R(k-1, j-1)] / (4^j - 1)
    static IntegrationResult romberg(
        core::RealFunction f,
        double a,
        double b,
        int maxOrder = 5,
        double tol = 1e-8,
        double trueVal = 0.0
    );

    // Gaussian Quadrature:
    // int_a^b f(x) dx = (b - a)/2 * sum_{i=1}^n w_i * f((b-a)/2 * t_i + (a+b)/2)
    // Supports orders 2, 3, 4, 5 with exact Gauss-Legendre nodes & weights
    static IntegrationResult gaussianQuadrature(
        core::RealFunction f,
        double a,
        double b,
        int numPoints = 3,
        double trueVal = 0.0
    );

    // Expression overloads
    static IntegrationResult integrateExpr(
        const std::string& expr,
        double a,
        double b,
        const std::string& method,
        int n = 100,
        double trueVal = 0.0
    );

    // Multi-method comparison table
    static std::string compareIntegrationMethods(
        core::RealFunction f,
        double a,
        double b,
        double trueVal
    );
};

} // namespace ees::unit2

