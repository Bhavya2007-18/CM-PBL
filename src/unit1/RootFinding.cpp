#include "ees/unit1/RootFinding.hpp"
#include "ees/core/ErrorAnalysis.hpp"
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>

namespace ees::unit1 {

core::NumericalResult RootFindingLab::bisection(
    core::RealFunction f,
    double a,
    double b,
    double tol,
    int maxIter
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Linear;
    result.convergenceRateExplanation = "Bisection halves the interval at each step: e_{k+1} = 0.5 * e_k (Order 1, Asymptotic error constant 0.5).";

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0) {
        result.success = false;
        result.status = core::ConvergenceStatus::InvalidInput;
        result.message = "Bracketing condition violated: f(a) and f(b) must have opposite signs. f(" +
                         std::to_string(a) + ")=" + std::to_string(fa) + ", f(" +
                         std::to_string(b) + ")=" + std::to_string(fb);
        return result;
    }

    if (std::abs(fa) < 1e-15) {
        result.success = true;
        result.solution = a;
        result.status = core::ConvergenceStatus::Converged;
        result.iterations = 0;
        result.message = "Exact root found at lower bracket a.";
        return result;
    }
    if (std::abs(fb) < 1e-15) {
        result.success = true;
        result.solution = b;
        result.status = core::ConvergenceStatus::Converged;
        result.iterations = 0;
        result.message = "Exact root found at upper bracket b.";
        return result;
    }

    double c = a;
    double prevC = a;

    for (int iter = 1; iter <= maxIter; ++iter) {
        c = 0.5 * (a + b);
        double fc = f(c);
        double err = (iter == 1) ? 0.5 * std::abs(b - a) : std::abs(c - prevC);
        double relErr = core::ErrorAnalysis::relativeError(c, prevC);

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.x = c;
        rec.fx = fc;
        rec.absoluteError = err;
        rec.relativeError = relErr;
        rec.aux["a"] = a;
        rec.aux["b"] = b;
        rec.aux["f(a)"] = fa;
        rec.aux["f(b)"] = fb;
        rec.aux["interval_width"] = std::abs(b - a);
        result.addRecord(rec);

        if (std::abs(fc) < 1e-15 || err < tol || 0.5 * std::abs(b - a) < tol) {
            result.success = true;
            result.solution = c;
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = err;
            result.finalRelativeError = relErr;
            result.message = "Bisection converged within tolerance after " + std::to_string(iter) + " iterations.";
            break;
        }

        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        prevC = c;
    }

    if (!result.success) {
        result.solution = c;
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = 0.5 * std::abs(b - a);
        result.message = "Maximum iterations reached before reaching target tolerance.";
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult RootFindingLab::newtonRaphson(
    core::RealFunction f,
    double x0,
    double tol,
    int maxIter,
    core::RealFunction df
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Quadratic;
    result.convergenceRateExplanation = "Newton-Raphson doubles the number of correct significant digits per iteration near a simple root (Order 2: e_{k+1} ≈ (f''(r)/(2f'(r))) * e_k^2).";

    double x = x0;

    for (int iter = 1; iter <= maxIter; ++iter) {
        double fx = f(x);
        double dfx = df ? df(x) : core::NumericalDifferentiation::centralDiff(f, x, 1e-6);

        if (std::abs(dfx) < 1e-15) {
            result.success = false;
            result.solution = x;
            result.status = core::ConvergenceStatus::DivisionByZero;
            result.message = "Derivative vanishes (f'(x) ≈ 0) at iteration " + std::to_string(iter) + ", causing division by zero.";
            break;
        }

        double xNext = x - (fx / dfx);
        double err = std::abs(xNext - x);
        double relErr = core::ErrorAnalysis::relativeError(xNext, x);

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.x = xNext;
        rec.fx = f(xNext);
        rec.absoluteError = err;
        rec.relativeError = relErr;
        rec.aux["f'(x)"] = dfx;
        rec.aux["step"] = -(fx / dfx);
        result.addRecord(rec);

        if (err < tol || std::abs(rec.fx) < tol) {
            result.success = true;
            result.solution = xNext;
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = err;
            result.finalRelativeError = relErr;
            result.message = "Newton-Raphson converged quadratically after " + std::to_string(iter) + " iterations.";
            break;
        }

        // Divergence detection
        if (!std::isfinite(xNext) || std::abs(xNext) > 1e12) {
            result.success = false;
            result.solution = xNext;
            result.status = core::ConvergenceStatus::Diverged;
            result.message = "Method diverged: iterates exceeded numeric stability threshold.";
            break;
        }

        x = xNext;
    }

    if (!result.success && result.status != core::ConvergenceStatus::DivisionByZero && result.status != core::ConvergenceStatus::Diverged) {
        result.solution = x;
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = result.history.empty() ? 0.0 : result.history.back().absoluteError;
        result.message = "Maximum iterations reached before reaching target tolerance.";
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult RootFindingLab::secant(
    core::RealFunction f,
    double x0,
    double x1,
    double tol,
    int maxIter
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Superlinear;
    result.convergenceRateExplanation = "Secant approximates the derivative via finite difference: order ≈ (1 + sqrt(5))/2 ≈ 1.618 (Superlinear).";

    double xPrev = x0;
    double xCurr = x1;
    double fPrev = f(xPrev);
    double fCurr = f(xCurr);

    for (int iter = 1; iter <= maxIter; ++iter) {
        double denom = fCurr - fPrev;
        if (std::abs(denom) < 1e-15) {
            result.success = false;
            result.solution = xCurr;
            result.status = core::ConvergenceStatus::DivisionByZero;
            result.message = "Zero secant slope: f(x_k) - f(x_{k-1}) ≈ 0 at iteration " + std::to_string(iter);
            break;
        }

        double xNext = xCurr - fCurr * (xCurr - xPrev) / denom;
        double fNext = f(xNext);
        double err = std::abs(xNext - xCurr);
        double relErr = core::ErrorAnalysis::relativeError(xNext, xCurr);

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.x = xNext;
        rec.fx = fNext;
        rec.absoluteError = err;
        rec.relativeError = relErr;
        rec.aux["x_prev"] = xPrev;
        rec.aux["slope"] = denom / (xCurr - xPrev);
        result.addRecord(rec);

        if (err < tol || std::abs(fNext) < tol) {
            result.success = true;
            result.solution = xNext;
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = err;
            result.finalRelativeError = relErr;
            result.message = "Secant method converged superlinearly after " + std::to_string(iter) + " iterations.";
            break;
        }

        if (!std::isfinite(xNext) || std::abs(xNext) > 1e12) {
            result.success = false;
            result.solution = xNext;
            result.status = core::ConvergenceStatus::Diverged;
            result.message = "Method diverged: iterates exceeded numeric stability threshold.";
            break;
        }

        xPrev = xCurr;
        fPrev = fCurr;
        xCurr = xNext;
        fCurr = fNext;
    }

    if (!result.success && result.status != core::ConvergenceStatus::DivisionByZero && result.status != core::ConvergenceStatus::Diverged) {
        result.solution = xCurr;
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = result.history.empty() ? 0.0 : result.history.back().absoluteError;
        result.message = "Maximum iterations reached before reaching target tolerance.";
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult RootFindingLab::bisectionExpr(
    const std::string& expr,
    double a,
    double b,
    double tol,
    int maxIter
) {
    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };
    return bisection(f, a, b, tol, maxIter);
}

core::NumericalResult RootFindingLab::newtonRaphsonExpr(
    const std::string& expr,
    double x0,
    double tol,
    int maxIter
) {
    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };
    return newtonRaphson(f, x0, tol, maxIter);
}

core::NumericalResult RootFindingLab::secantExpr(
    const std::string& expr,
    double x0,
    double x1,
    double tol,
    int maxIter
) {
    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };
    return secant(f, x0, x1, tol, maxIter);
}

std::string RootFindingLab::compareRootMethods(
    core::RealFunction f,
    double a,
    double b,
    double x0,
    double tol,
    int maxIter
) {
    auto resBis = bisection(f, a, b, tol, maxIter);
    auto resNew = newtonRaphson(f, x0, tol, maxIter);
    auto resSec = secant(f, a, b, tol, maxIter);

    std::ostringstream ss;
    ss << "========================================================================================\n";
    ss << "                            ROOT-FINDING METHOD COMPARISON                              \n";
    ss << "========================================================================================\n";
    ss << std::left << std::setw(16) << "Method"
       << std::setw(16) << "Computed Root"
       << std::setw(12) << "Iterations"
       << std::setw(16) << "Final Error"
       << std::setw(14) << "Time (ms)"
       << std::setw(14) << "Convergence" << "\n";
    ss << std::string(88, '-') << "\n";

    auto printRow = [&ss](const std::string& name, const core::NumericalResult& r) {
        ss << std::left << std::setw(16) << name
           << std::setprecision(8) << std::fixed << std::setw(16) << r.solution
           << std::setw(12) << r.iterations
           << std::scientific << std::setprecision(4) << std::setw(16) << r.finalError
           << std::fixed << std::setprecision(3) << std::setw(14) << r.executionTimeMs
           << std::setw(14) << core::convergenceRateToString(r.convergenceRate) << "\n";
    };

    printRow("Bisection", resBis);
    printRow("Newton-Raphson", resNew);
    printRow("Secant", resSec);
    ss << std::string(88, '-') << "\n";

    ss << "ANALYSIS & TRADE-OFFS:\n";
    ss << "- Bisection: Guaranteed to converge if bracket f(a)*f(b) < 0 exists, but requires the most iterations (linear).\n";
    ss << "- Newton-Raphson: Fastest local convergence (quadratic, ~doubles digits/iter), requires continuous derivative.\n";
    ss << "- Secant: Almost as fast as Newton (order ~1.618), derivative-free, slightly more sensitive to starting guesses.\n";

    return ss.str();
}

} // namespace ees::unit1

