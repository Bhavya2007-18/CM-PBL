#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/NumericalResult.hpp"
#include "ees/core/ExpressionParser.hpp"
#include <string>

namespace ees::unit1 {

struct RootFindingParams {
    double a{0.0};              // Lower bracket or initial guess 1
    double b{0.0};              // Upper bracket or initial guess 2
    double x0{0.0};             // Single starting guess (for Newton-Raphson)
    double tolerance{1e-7};     // Target convergence tolerance
    int maxIterations{100};     // Max iteration guard
    bool verbose{false};
};

class RootFindingLab {
public:
    // Bisection Method:
    // Requires f(a) * f(b) < 0. Convergence is guaranteed and linear.
    static core::NumericalResult bisection(
        core::RealFunction f,
        double a,
        double b,
        double tol = 1e-7,
        int maxIter = 100
    );

    // Newton-Raphson Method:
    // x_{n+1} = x_n - f(x_n)/f'(x_n). Quadratic convergence near simple roots.
    static core::NumericalResult newtonRaphson(
        core::RealFunction f,
        double x0,
        double tol = 1e-7,
        int maxIter = 100,
        core::RealFunction df = nullptr
    );

    // Secant Method:
    // x_{n+1} = x_n - f(x_n)*(x_n - x_{n-1}) / (f(x_n) - f(x_{n-1})).
    // Superlinear convergence (order ~1.618), derivative-free.
    static core::NumericalResult secant(
        core::RealFunction f,
        double x0,
        double x1,
        double tol = 1e-7,
        int maxIter = 100
    );

    // Expression-based overloads
    static core::NumericalResult bisectionExpr(
        const std::string& expr,
        double a,
        double b,
        double tol = 1e-7,
        int maxIter = 100
    );

    static core::NumericalResult newtonRaphsonExpr(
        const std::string& expr,
        double x0,
        double tol = 1e-7,
        int maxIter = 100
    );

    static core::NumericalResult secantExpr(
        const std::string& expr,
        double x0,
        double x1,
        double tol = 1e-7,
        int maxIter = 100
    );

    // Method comparison on the same problem
    static std::string compareRootMethods(
        core::RealFunction f,
        double a,
        double b,
        double x0,
        double tol = 1e-7,
        int maxIter = 100
    );
};

} // namespace ees::unit1

