#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include "ees/core/NumericalResult.hpp"
#include "ees/core/ExpressionParser.hpp"
#include <vector>
#include <functional>

namespace ees::unit1 {

class OptimizationLab {
public:
    // 1D Optimization: Fibonacci Search
    // Given unimodal function f on [a, b], shrinks interval using Fibonacci sequence
    static core::NumericalResult fibonacciSearch(
        core::RealFunction f,
        double a,
        double b,
        double tol = 1e-5,
        int maxIter = 40,
        core::OptimizationGoal goal = core::OptimizationGoal::Minimize
    );

    // 1D Optimization: Golden Section Search
    // Uses golden ratio phi = (sqrt(5) - 1)/2 ≈ 0.6180339887
    static core::NumericalResult goldenSectionSearch(
        core::RealFunction f,
        double a,
        double b,
        double tol = 1e-5,
        int maxIter = 100,
        core::OptimizationGoal goal = core::OptimizationGoal::Minimize
    );

    // 1D Optimization: Newton's Method
    // Finds stationary point: x_{k+1} = x_k - f'(x_k) / f''(x_k)
    static core::NumericalResult newtonsMethod1D(
        core::RealFunction f,
        double x0,
        double tol = 1e-7,
        int maxIter = 100
    );

    // Multivariate Optimization: Steepest Descent (Gradient Descent)
    // x_{k+1} = x_k - alpha * grad(f(x_k)) with backtracking line search
    static core::NumericalResult steepestDescent(
        const std::function<double(const core::Vector&)>& f,
        const core::Vector& initialPoint,
        double stepSize = 0.1,
        double tol = 1e-5,
        int maxIter = 200,
        bool useBacktracking = true
    );

    // Multivariate Optimization: Nelder-Mead Simplex Algorithm
    // Direct search method based on simplex transformations:
    // Reflection -> Expansion / Contraction -> Shrink
    static core::NumericalResult nelderMead(
        const std::function<double(const core::Vector&)>& f,
        const core::Vector& initialPoint,
        double initialStep = 0.5,
        double tol = 1e-5,
        int maxIter = 300
    );
};

} // namespace ees::unit1

