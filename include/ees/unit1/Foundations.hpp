#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/NumericalResult.hpp"
#include "ees/core/ExpressionParser.hpp"
#include <vector>
#include <string>

namespace ees::unit1 {

struct TaylorTerm {
    int order{0};
    double derivativeVal{0.0};
    double termVal{0.0};
    double partialSum{0.0};
    double absError{0.0};
    double relError{0.0};
    double pctError{0.0};
};

struct TaylorResult {
    double exactValue{0.0};
    double approxValue{0.0};
    double x0{0.0};
    double x{0.0};
    int numTerms{0};
    std::vector<TaylorTerm> terms;
    std::string summary() const;
};

struct TheoremResult {
    bool theoremHolds{false};
    std::string theoremName;
    double a{0.0};
    double b{0.0};
    double fa{0.0};
    double fb{0.0};
    double targetSlope{0.0}; // 0 for Rolle's, [f(b)-f(a)]/(b-a) for MVT
    std::vector<double> cValues; // Points in (a, b) where f'(c) == targetSlope
    std::string explanation;
};

class FoundationsLab {
public:
    // Taylor Series expansion of f(x) about x0 evaluated at x up to n terms
    static TaylorResult computeTaylorSeries(
        core::RealFunction f,
        double x0,
        double x,
        int maxOrder,
        double trueValue
    );

    static TaylorResult computeTaylorSeriesExpr(
        const std::string& expression,
        double x0,
        double x,
        int maxOrder,
        double trueValue
    );

    // Rolle's Theorem: If f continuous on [a,b], diff on (a,b), and f(a)==f(b), then exists c in (a,b) where f'(c)=0
    static TheoremResult verifyRollesTheorem(
        core::RealFunction f,
        double a,
        double b,
        double tol = 1e-5
    );

    // Mean Value Theorem: If f continuous on [a,b], diff on (a,b), then exists c in (a,b) where f'(c) = (f(b)-f(a))/(b-a)
    static TheoremResult verifyMeanValueTheorem(
        core::RealFunction f,
        double a,
        double b,
        double tol = 1e-5
    );

    // Computer Arithmetic / Loss of Significance Lab report
    static std::string runArithmeticStabilityLab();
};

} // namespace ees::unit1

