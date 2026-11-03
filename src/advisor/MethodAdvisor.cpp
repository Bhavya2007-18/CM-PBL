#include "ees/advisor/MethodAdvisor.hpp"
#include <sstream>

namespace ees::advisor {

std::string Recommendation::toString() const {
    std::ostringstream ss;
    ss << "=== Numerical Method Advisor Guidance ===\n";
    ss << "Problem Characteristics Detected:\n";
    for (const auto& ch : detectedCharacteristics) {
        ss << "  [✓] " << ch << "\n";
    }
    ss << "\nRECOMMENDED PRIMARY METHOD: " << primaryMethod << "\n";
    if (!alternateMethod.empty()) {
        ss << "BACKUP / ALTERNATIVE METHOD: " << alternateMethod << "\n";
    }
    ss << "\nWHY THIS METHOD WAS RECOMMENDED:\n" << mathematicalJustification << "\n";
    if (!cautionaryNotes.empty()) {
        ss << "\nCAUTIONS & LIMITATIONS:\n" << cautionaryNotes << "\n";
    }
    return ss.str();
}

Recommendation MethodAdvisor::advise(const ProblemProfile& profile) {
    Recommendation rec;

    switch (profile.category) {
        case ProblemCategory::RootFinding:
            if (profile.hasBracketingInterval) {
                rec.detectedCharacteristics.push_back("Bracketing interval [a, b] with opposite signs available");
            }
            if (profile.hasAnalyticDerivative) {
                rec.detectedCharacteristics.push_back("Continuous derivative f'(x) readily computable");
            } else {
                rec.detectedCharacteristics.push_back("Derivative f'(x) absent or computationally prohibitive");
            }

            if (profile.hasAnalyticDerivative) {
                rec.primaryMethod = "Newton-Raphson Method";
                rec.alternateMethod = "Secant Method";
                rec.mathematicalJustification =
                    "Newton-Raphson yields asymptotic quadratic convergence (Order 2), doubling the number of correct\n"
                    "significant digits per cycle near simple roots.";
                rec.cautionaryNotes =
                    "Requires starting guess sufficiently close to root; fails if f'(x) approaches zero (inflection point).";
            } else if (profile.hasBracketingInterval) {
                rec.primaryMethod = "Bisection Method";
                rec.alternateMethod = "Secant Method";
                rec.mathematicalJustification =
                    "Bisection guarantees convergence via the Intermediate Value Theorem regardless of function shape,\n"
                    "halving the uncertainty interval at every iteration.";
                rec.cautionaryNotes =
                    "Convergence is linear (Order 1); requires more iterations than open methods.";
            } else {
                rec.primaryMethod = "Secant Method";
                rec.alternateMethod = "Bisection (upon locating sign change)";
                rec.mathematicalJustification =
                    "Secant achieves superlinear convergence (Order ~1.618) without requiring explicit derivative functions.";
                rec.cautionaryNotes =
                    "Can diverge if two consecutive approximations yield equal function values.";
            }
            break;

        case ProblemCategory::LinearSystem:
            if (profile.isMatrixSymmetricPositiveDefinite) {
                rec.detectedCharacteristics.push_back("Matrix A is Symmetric and Positive-Definite (x^T A x > 0)");
                rec.primaryMethod = "Cholesky Factorization (A = L * L^T)";
                rec.alternateMethod = "Doolittle LU Factorization";
                rec.mathematicalJustification =
                    "Cholesky requires roughly half the computational effort of standard LU decomposition (n^3 / 3 FLOPs vs 2n^3 / 3)\n"
                    "and is unconditionally numerically stable without row pivoting.";
            } else {
                rec.detectedCharacteristics.push_back("General square nonsymmetric matrix");
                rec.primaryMethod = "Gauss Elimination with Partial Pivoting";
                rec.alternateMethod = "Doolittle LU Decomposition";
                rec.mathematicalJustification =
                    "Partial pivoting bounds multiplier magnitudes (|m_{ik}| <= 1), suppressing round-off magnification.";
                rec.cautionaryNotes =
                    "Ensure matrix determinant is non-zero to avoid singular system exceptions.";
            }
            break;

        case ProblemCategory::NumericalIntegration:
            rec.detectedCharacteristics.push_back("Smooth integrand over finite limits [a, b]");
            if (profile.isSmoothContinuous) {
                rec.primaryMethod = "Gauss-Legendre Quadrature";
                rec.alternateMethod = "Romberg Extrapolation";
                rec.mathematicalJustification =
                    "An n-point Gauss-Legendre quadrature exactly integrates any polynomial of degree up to 2n - 1,\n"
                    "delivering maximum accuracy per function evaluation.";
            } else {
                rec.primaryMethod = "Composite Simpson's 1/3 Rule";
                rec.alternateMethod = "Trapezoidal Rule";
                rec.mathematicalJustification =
                    "Simpson's 1/3 fits parabolic arcs across paired subintervals with global truncation error O(h^4).";
            }
            break;

        case ProblemCategory::Interpolation:
            if (profile.isEquispacedData) {
                rec.detectedCharacteristics.push_back("Equispaced abscissas (constant grid step h)");
                rec.primaryMethod = "Gregory-Newton Forward / Backward Interpolation";
                rec.alternateMethod = "Natural Cubic Spline";
                rec.mathematicalJustification =
                    "Finite differences simplify polynomial evaluation and provide direct inspection of higher-order derivative trends.";
            } else {
                rec.detectedCharacteristics.push_back("Arbitrary / unequally spaced data points");
                rec.primaryMethod = "Newton Divided Difference / Natural Cubic Spline";
                rec.alternateMethod = "Lagrange Polynomial";
                rec.mathematicalJustification =
                    "Divided differences allow incremental addition of data knots without recomputing existing basis terms.";
            }
            break;

        case ProblemCategory::ODEInitialValue:
            rec.detectedCharacteristics.push_back("First-order Initial Value Problem dy/dx = f(x, y)");
            rec.primaryMethod = "Classic Runge-Kutta 4th-Order (RK4)";
            rec.alternateMethod = "Adams-Bashforth-Moulton Predictor-Corrector";
            rec.mathematicalJustification =
                "RK4 is self-starting, achieves 4th-order global accuracy O(h^4), and possesses a broad linear stability domain.";
            break;

        default:
            rec.primaryMethod = "Standard Syllabus Algorithm";
            rec.mathematicalJustification = "Default heuristic based on standard computational methods.";
            break;
    }

    return rec;
}

} // namespace ees::advisor

