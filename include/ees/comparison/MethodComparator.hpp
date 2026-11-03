#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/NumericalResult.hpp"
#include <string>
#include <vector>

namespace ees::comparison {

struct MethodBenchmarkRow {
    std::string methodName;
    double numericalResult{0.0};
    int iterations{0};
    double finalError{0.0};
    double executionTimeMs{0.0};
    std::string convergenceType;
    std::string stabilityRemarks;
};

struct BenchmarkReport {
    std::string problemDescription;
    std::vector<MethodBenchmarkRow> rows;
    std::string recommendedMethod;
    std::string tradeOffAnalysis;
    std::string summary() const;
};

class MethodComparator {
public:
    // Compare Root Finding Methods on f(x) = 0
    static BenchmarkReport compareRoots(
        core::RealFunction f,
        double a,
        double b,
        double x0,
        double tol = 1e-7,
        const std::string& eqStr = "f(x) = 0"
    );

    // Compare 1D Optimization Methods on f(x)
    static BenchmarkReport compare1DOptimization(
        core::RealFunction f,
        double a,
        double b,
        double x0,
        double tol = 1e-5,
        const std::string& funcStr = "f(x)"
    );

    // Compare Quadrature Methods on int_a^b f(x) dx
    static BenchmarkReport compareIntegration(
        core::RealFunction f,
        double a,
        double b,
        double trueValue,
        const std::string& funcStr = "f(x)"
    );

    // Compare ODE Solvers on dy/dx = f(x, y)
    static BenchmarkReport compareODESolvers(
        std::function<double(double, double)> f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol,
        const std::string& odeStr = "dy/dx = f(x, y)"
    );
};

} // namespace ees::comparison

