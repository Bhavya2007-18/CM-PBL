#include "ees/unit1/Foundations.hpp"
#include "ees/core/ErrorAnalysis.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <functional>

namespace ees::unit1 {

std::string TaylorResult::summary() const {
    std::ostringstream ss;
    ss << "=== Taylor Series Approximation ===\n";
    ss << "Expansion Point (x0): " << x0 << "\n";
    ss << "Evaluation Point (x):  " << x << "\n";
    ss << "True Value:            " << std::setprecision(10) << exactValue << "\n";
    ss << "Taylor Approx (Order " << numTerms - 1 << "): " << approxValue << "\n";
    ss << "Final Absolute Error:  " << std::scientific << std::abs(exactValue - approxValue) << "\n";
    ss << "Final Relative Error:  " << core::ErrorAnalysis::relativeError(exactValue, approxValue) << "\n\n";

    ss << std::fixed << std::setprecision(6);
    ss << std::setw(6) << "Order" << " | "
       << std::setw(14) << "Term Value" << " | "
       << std::setw(14) << "Partial Sum" << " | "
       << std::setw(14) << "Abs Error" << " | "
       << std::setw(12) << "% Error" << "\n";
    ss << std::string(70, '-') << "\n";

    for (const auto& t : terms) {
        ss << std::setw(6) << t.order << " | "
           << std::setw(14) << t.termVal << " | "
           << std::setw(14) << t.partialSum << " | "
           << std::setw(14) << t.absError << " | "
           << std::setw(12) << t.pctError << "\n";
    }
    return ss.str();
}

TaylorResult FoundationsLab::computeTaylorSeries(
    core::RealFunction f,
    double x0,
    double x,
    int maxOrder,
    double trueValue
) {
    TaylorResult res;
    res.x0 = x0;
    res.x = x;
    res.numTerms = maxOrder + 1;
    res.exactValue = trueValue;

    double h = x - x0;
    double currentSum = 0.0;
    double factorial = 1.0;
    double hPow = 1.0;

    // Helper for higher-order derivatives using finite differences
    auto computeNthDerivative = [&](int n, double pt) -> double {
        if (n == 0) return f(pt);
        if (n == 1) return core::NumericalDifferentiation::centralDiff(f, pt, 1e-5);
        if (n == 2) return core::NumericalDifferentiation::secondDiff(f, pt, 1e-4);

        // General finite difference stencil for higher order derivatives:
        // sum_{k=0}^n (-1)^(n-k) * binom(n, k) * f(x + k*step) / step^n
        double step = 0.05;
        double sum = 0.0;
        std::vector<double> binom(n + 1, 1.0);
        for (int i = 1; i <= n; ++i) {
            for (int j = i - 1; j > 0; --j) {
                binom[j] += binom[j - 1];
            }
        }
        for (int k = 0; k <= n; ++k) {
            double sign = ((n - k) % 2 == 0) ? 1.0 : -1.0;
            sum += sign * binom[k] * f(pt + (k - n * 0.5) * step);
        }
        return sum / std::pow(step, n);
    };

    for (int n = 0; n <= maxOrder; ++n) {
        if (n > 0) {
            factorial *= n;
            hPow *= h;
        }
        double dVal = computeNthDerivative(n, x0);
        double termVal = (dVal / factorial) * hPow;
        currentSum += termVal;

        TaylorTerm term;
        term.order = n;
        term.derivativeVal = dVal;
        term.termVal = termVal;
        term.partialSum = currentSum;
        term.absError = core::ErrorAnalysis::absoluteError(trueValue, currentSum);
        term.relError = core::ErrorAnalysis::relativeError(trueValue, currentSum);
        term.pctError = core::ErrorAnalysis::percentageError(trueValue, currentSum);
        res.terms.push_back(term);
    }

    res.approxValue = currentSum;
    return res;
}

TaylorResult FoundationsLab::computeTaylorSeriesExpr(
    const std::string& expression,
    double x0,
    double x,
    int maxOrder,
    double trueValue
) {
    core::ExpressionParser parser(expression);
    auto f = [parser](double val) mutable {
        return parser.evaluate(val);
    };
    return computeTaylorSeries(f, x0, x, maxOrder, trueValue);
}

TheoremResult FoundationsLab::verifyRollesTheorem(
    core::RealFunction f,
    double a,
    double b,
    double tol
) {
    TheoremResult res;
    res.theoremName = "Rolle's Theorem";
    res.a = a;
    res.b = b;
    res.fa = f(a);
    res.fb = f(b);
    res.targetSlope = 0.0;

    std::ostringstream ss;
    ss << "Verification of Rolle's Theorem on [" << a << ", " << b << "]:\n";
    ss << "1. Hypothesis f(a) == f(b): f(" << a << ") = " << res.fa << ", f(" << b << ") = " << res.fb << "\n";

    if (std::abs(res.fa - res.fb) > 1e-4) {
        res.theoremHolds = false;
        ss << "   -> Rolle's Theorem preconditions NOT met: f(a) != f(b).\n";
        res.explanation = ss.str();
        return res;
    }

    ss << "   -> Preconditions satisfied: f(a) == f(b).\n";
    ss << "2. Searching for critical point c in (" << a << ", " << b << ") where f'(c) = 0...\n";

    // Subdivide interval to locate sign changes of derivative f'(x)
    const int numDivs = 1000;
    double dx = (b - a) / numDivs;
    auto df = [f](double x) {
        return core::NumericalDifferentiation::centralDiff(f, x, 1e-5);
    };

    double prevX = a + dx;
    double prevDf = df(prevX);

    for (int i = 1; i < numDivs - 1; ++i) {
        double curX = a + i * dx;
        double curDf = df(curX);

        if (std::abs(curDf) < tol) {
            res.cValues.push_back(curX);
        } else if (prevDf * curDf < 0.0) {
            // Sign change bracket for root of f'(x)
            double low = prevX, high = curX;
            for (int it = 0; it < 30; ++it) {
                double mid = 0.5 * (low + high);
                double dMid = df(mid);
                if (std::abs(dMid) < tol) {
                    low = mid;
                    break;
                }
                if (df(low) * dMid < 0.0) high = mid;
                else low = mid;
            }
            res.cValues.push_back(low);
        }
        prevX = curX;
        prevDf = curDf;
    }

    if (!res.cValues.empty()) {
        res.theoremHolds = true;
        ss << "   -> Found " << res.cValues.size() << " point(s) c in (" << a << ", " << b << ") with f'(c) ≈ 0:\n";
        for (double c : res.cValues) {
            ss << "      c = " << c << " (f'(c) = " << df(c) << ")\n";
        }
    } else {
        res.theoremHolds = false;
        ss << "   -> Warning: No stationary point located in interior (function may not be continuously differentiable).\n";
    }

    res.explanation = ss.str();
    return res;
}

TheoremResult FoundationsLab::verifyMeanValueTheorem(
    core::RealFunction f,
    double a,
    double b,
    double tol
) {
    TheoremResult res;
    res.theoremName = "Mean Value Theorem (Lagrange)";
    res.a = a;
    res.b = b;
    res.fa = f(a);
    res.fb = f(b);
    res.targetSlope = (res.fb - res.fa) / (b - a);

    std::ostringstream ss;
    ss << "Verification of Mean Value Theorem on [" << a << ", " << b << "]:\n";
    ss << "Secant slope [f(b) - f(a)] / (b - a) = [" << res.fb << " - " << res.fa << "] / " << (b - a)
       << " = " << res.targetSlope << "\n";
    ss << "Searching for point(s) c in (" << a << ", " << b << ") such that f'(c) = " << res.targetSlope << "...\n";

    auto g = [f, res](double x) {
        return core::NumericalDifferentiation::centralDiff(f, x, 1e-5) - res.targetSlope;
    };

    const int numDivs = 1000;
    double dx = (b - a) / numDivs;
    double prevX = a + dx;
    double prevG = g(prevX);

    for (int i = 1; i < numDivs - 1; ++i) {
        double curX = a + i * dx;
        double curG = g(curX);

        if (std::abs(curG) < tol) {
            res.cValues.push_back(curX);
        } else if (prevG * curG < 0.0) {
            double low = prevX, high = curX;
            for (int it = 0; it < 30; ++it) {
                double mid = 0.5 * (low + high);
                double gMid = g(mid);
                if (std::abs(gMid) < tol) {
                    low = mid;
                    break;
                }
                if (g(low) * gMid < 0.0) high = mid;
                else low = mid;
            }
            res.cValues.push_back(low);
        }
        prevX = curX;
        prevG = curG;
    }

    if (!res.cValues.empty()) {
        res.theoremHolds = true;
        ss << "-> Found " << res.cValues.size() << " point(s) c in (" << a << ", " << b << ") where instantaneous tangent == secant slope:\n";
        for (double c : res.cValues) {
            double actualSlope = core::NumericalDifferentiation::centralDiff(f, c, 1e-5);
            ss << "   c = " << c << " (f'(c) = " << actualSlope << ", target = " << res.targetSlope << ")\n";
        }
    } else {
        res.theoremHolds = false;
        ss << "-> Notice: No point located in interior.\n";
    }

    res.explanation = ss.str();
    return res;
}

std::string FoundationsLab::runArithmeticStabilityLab() {
    std::ostringstream ss;
    ss << "========================================================================\n";
    ss << "       COMPUTER ARITHMETIC & NUMERICAL STABILITY LABORATORY             \n";
    ss << "========================================================================\n\n";

    ss << "1. CATASTROPHIC CANCELLATION DEMONSTRATION:\n";
    ss << "   Problem: Evaluate f(x) = sqrt(x + 1) - sqrt(x) for large x.\n";
    ss << "   Mathematically, as x -> inf, subtracting nearly equal floating-point\n";
    ss << "   numbers cancels significant higher mantissa bits, causing severe round-off error.\n";
    ss << "   Algebraically equivalent stable form: 1.0 / (sqrt(x + 1) + sqrt(x))\n\n";

    double testVals[] = {1.0, 100.0, 1e6, 1e10, 1e14, 1e16};
    ss << std::setw(12) << "x" << " | "
       << std::setw(20) << "Direct f(x)" << " | "
       << std::setw(20) << "Stable f(x)" << " | "
       << std::setw(18) << "Lost Sig Digits" << "\n";
    ss << std::string(78, '-') << "\n";

    for (double x : testVals) {
        double d = core::ErrorAnalysis::cancellationDirect(x);
        double s = core::ErrorAnalysis::cancellationStable(x);
        double relErr = core::ErrorAnalysis::relativeError(s, d);
        int sigDigits = core::ErrorAnalysis::significantDigits(relErr);

        ss << std::scientific << std::setprecision(4)
           << std::setw(12) << x << " | "
           << std::setprecision(10) << std::setw(20) << d << " | "
           << std::setw(20) << s << " | "
           << std::setw(18) << (16 - sigDigits) << "\n";
    }

    ss << "\n2. TRIGONOMETRIC LOSS OF SIGNIFICANCE FOR SMALL x:\n";
    ss << "   Problem: Evaluate g(x) = (1 - cos(x)) / x^2 as x -> 0.\n";
    ss << "   For very small x, cos(x) rounds to 1.0 in IEEE-754, causing (1 - 1)/x^2 = 0!\n";
    ss << "   Stable identity: 2*sin^2(x/2)/x^2 = 0.5 * (sin(x/2)/(x/2))^2\n\n";

    double smallVals[] = {1.0, 0.1, 1e-4, 1e-7, 1e-8, 1e-9};
    ss << std::setw(12) << "x" << " | "
       << std::setw(20) << "Direct g(x)" << " | "
       << std::setw(20) << "Stable g(x)" << " | "
       << std::setw(18) << "Theoretical" << "\n";
    ss << std::string(78, '-') << "\n";

    for (double x : smallVals) {
        double gd = core::ErrorAnalysis::cancellationTrigDirect(x);
        double gs = core::ErrorAnalysis::cancellationTrigStable(x);
        ss << std::scientific << std::setprecision(4)
           << std::setw(12) << x << " | "
           << std::setprecision(10) << std::setw(20) << gd << " | "
           << std::setw(20) << gs << " | "
           << std::setw(18) << 0.5 << "\n";
    }

    ss << "\nConclusion: In scientific software, mathematically identical formulas\n";
    ss << "produce dramatically distinct results depending on floating-point conditioning.\n";
    return ss.str();
}

} // namespace ees::unit1

