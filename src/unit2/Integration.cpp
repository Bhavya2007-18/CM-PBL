#include "ees/unit2/Integration.hpp"
#include "ees/core/ErrorAnalysis.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace ees::unit2 {

std::string IntegrationResult::summary() const {
    std::ostringstream ss;
    ss << "=== Numerical Quadrature Result (" << method << ") ===\n";
    ss << "Computed Integral:  " << std::setprecision(10) << value << "\n";
    if (std::abs(trueValue) > 1e-15) {
        ss << "True Value:         " << trueValue << "\n";
        ss << "Absolute Error:     " << std::scientific << absError << "\n";
        ss << "Relative Error:     " << relError << "\n";
    }
    ss << "Subintervals:       " << subintervals << "\n";
    ss << "Function Calls:     " << functionEvaluations << "\n";

    if (!rombergTable.empty()) {
        ss << "\nRomberg Table (Extrapolations):\n";
        for (size_t r = 0; r < rombergTable.rows(); ++r) {
            ss << "R(" << r << "): ";
            for (size_t c = 0; c <= r; ++c) {
                ss << std::fixed << std::setprecision(8) << rombergTable(r, c) << "  ";
            }
            ss << "\n";
        }
    }
    return ss.str();
}

IntegrationResult IntegrationLab::trapezoidal(
    core::RealFunction f,
    double a,
    double b,
    int n,
    double trueVal
) {
    if (n < 1) throw std::invalid_argument("Subintervals must be >= 1 for Trapezoidal rule");

    double h = (b - a) / static_cast<double>(n);
    double sum = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    double val = sum * h;

    IntegrationResult res;
    res.value = val;
    res.trueValue = trueVal;
    res.subintervals = n;
    res.functionEvaluations = n + 1;
    res.method = "Trapezoidal Rule";
    res.absError = core::ErrorAnalysis::absoluteError(trueVal, val);
    res.relError = core::ErrorAnalysis::relativeError(trueVal, val);
    return res;
}

IntegrationResult IntegrationLab::simpson13(
    core::RealFunction f,
    double a,
    double b,
    int n,
    double trueVal
) {
    if (n % 2 != 0) {
        // Enforce even number of subintervals for Simpson's 1/3 rule
        n += 1;
    }

    double h = (b - a) / static_cast<double>(n);
    double sum = f(a) + f(b);

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2 == 1) {
            sum += 4.0 * f(x);
        } else {
            sum += 2.0 * f(x);
        }
    }
    double val = (h / 3.0) * sum;

    IntegrationResult res;
    res.value = val;
    res.trueValue = trueVal;
    res.subintervals = n;
    res.functionEvaluations = n + 1;
    res.method = "Simpson's 1/3 Rule";
    res.absError = core::ErrorAnalysis::absoluteError(trueVal, val);
    res.relError = core::ErrorAnalysis::relativeError(trueVal, val);
    return res;
}

IntegrationResult IntegrationLab::simpson38(
    core::RealFunction f,
    double a,
    double b,
    int n,
    double trueVal
) {
    // n must be a multiple of 3
    if (n % 3 != 0) {
        n += (3 - (n % 3));
    }

    double h = (b - a) / static_cast<double>(n);
    double sum = f(a) + f(b);

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 3 == 0) {
            sum += 2.0 * f(x);
        } else {
            sum += 3.0 * f(x);
        }
    }
    double val = (3.0 * h / 8.0) * sum;

    IntegrationResult res;
    res.value = val;
    res.trueValue = trueVal;
    res.subintervals = n;
    res.functionEvaluations = n + 1;
    res.method = "Simpson's 3/8 Rule";
    res.absError = core::ErrorAnalysis::absoluteError(trueVal, val);
    res.relError = core::ErrorAnalysis::relativeError(trueVal, val);
    return res;
}

IntegrationResult IntegrationLab::romberg(
    core::RealFunction f,
    double a,
    double b,
    int maxOrder,
    double tol,
    double trueVal
) {
    if (maxOrder < 2) maxOrder = 2;
    core::Matrix R = core::Matrix::zeros(maxOrder, maxOrder);

    int totalEvals = 0;

    // R(0, 0)
    double h = b - a;
    R(0, 0) = 0.5 * h * (f(a) + f(b));
    totalEvals += 2;

    double finalVal = R(0, 0);

    for (int k = 1; k < maxOrder; ++k) {
        // Compute trapezoidal rule with 2^k subintervals
        int n = 1 << k; // 2^k
        double hk = (b - a) / static_cast<double>(n);

        // Sum over newly added odd points
        double sumNew = 0.0;
        for (int i = 1; i < n; i += 2) {
            sumNew += f(a + i * hk);
            totalEvals++;
        }
        R(k, 0) = 0.5 * R(k - 1, 0) + hk * sumNew;

        // Richardson extrapolation column progression:
        // R(k, j) = R(k, j-1) + (R(k, j-1) - R(k-1, j-1)) / (4^j - 1)
        for (int j = 1; j <= k; ++j) {
            double factor = std::pow(4.0, j);
            R(k, j) = R(k, j - 1) + (R(k, j - 1) - R(k - 1, j - 1)) / (factor - 1.0);
        }

        finalVal = R(k, k);
        if (std::abs(R(k, k) - R(k - 1, k - 1)) < tol) {
            break;
        }
    }

    IntegrationResult res;
    res.value = finalVal;
    res.trueValue = trueVal;
    res.subintervals = 1 << (maxOrder - 1);
    res.functionEvaluations = totalEvals;
    res.method = "Romberg Extrapolation";
    res.rombergTable = R;
    res.absError = core::ErrorAnalysis::absoluteError(trueVal, finalVal);
    res.relError = core::ErrorAnalysis::relativeError(trueVal, finalVal);
    return res;
}

IntegrationResult IntegrationLab::gaussianQuadrature(
    core::RealFunction f,
    double a,
    double b,
    int numPoints,
    double trueVal
) {
    // Gauss-Legendre quadrature points (t_i in [-1, 1]) and weights (w_i)
    std::vector<double> t;
    std::vector<double> w;

    if (numPoints == 2) {
        double val = 1.0 / std::sqrt(3.0); // 0.577350269189626
        t = {-val, val};
        w = {1.0, 1.0};
    } else if (numPoints == 3) {
        double val = std::sqrt(3.0 / 5.0); // 0.774596669241483
        t = {-val, 0.0, val};
        w = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    } else if (numPoints == 4) {
        double t1 = std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
        double t2 = std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
        double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
        double w2 = (18.0 - std::sqrt(30.0)) / 36.0;
        t = {-t2, -t1, t1, t2};
        w = {w2, w1, w1, w2};
    } else if (numPoints == 5) {
        double t1 = (1.0 / 3.0) * std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0));
        double t2 = (1.0 / 3.0) * std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0));
        double w0 = 128.0 / 225.0;
        double w1 = (322.0 + 13.0 * std::sqrt(70.0)) / 900.0;
        double w2 = (322.0 - 13.0 * std::sqrt(70.0)) / 900.0;
        t = {-t2, -t1, 0.0, t1, t2};
        w = {w2, w1, w0, w1, w2};
    } else {
        throw std::invalid_argument("Gaussian Quadrature supports orders 2, 3, 4, and 5.");
    }

    // Affine variable substitution: x = (b - a)/2 * t + (a + b)/2, dx = (b - a)/2 * dt
    double halfLen = 0.5 * (b - a);
    double mid = 0.5 * (a + b);
    double sum = 0.0;

    for (size_t i = 0; i < t.size(); ++i) {
        double x = halfLen * t[i] + mid;
        sum += w[i] * f(x);
    }
    double val = halfLen * sum;

    IntegrationResult res;
    res.value = val;
    res.trueValue = trueVal;
    res.subintervals = numPoints;
    res.functionEvaluations = numPoints;
    res.method = "Gauss-Legendre (" + std::to_string(numPoints) + "-point)";
    res.absError = core::ErrorAnalysis::absoluteError(trueVal, val);
    res.relError = core::ErrorAnalysis::relativeError(trueVal, val);
    return res;
}

IntegrationResult IntegrationLab::integrateExpr(
    const std::string& expr,
    double a,
    double b,
    const std::string& method,
    int n,
    double trueVal
) {
    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };

    if (method == "trapezoidal") return trapezoidal(f, a, b, n, trueVal);
    if (method == "simpson13")   return simpson13(f, a, b, n, trueVal);
    if (method == "simpson38")   return simpson38(f, a, b, n, trueVal);
    if (method == "romberg")     return romberg(f, a, b, 5, 1e-8, trueVal);
    if (method == "gaussian")    return gaussianQuadrature(f, a, b, 3, trueVal);

    throw std::invalid_argument("Unknown integration method: " + method);
}

std::string IntegrationLab::compareIntegrationMethods(
    core::RealFunction f,
    double a,
    double b,
    double trueVal
) {
    auto rTrap = trapezoidal(f, a, b, 100, trueVal);
    auto rS13  = simpson13(f, a, b, 100, trueVal);
    auto rS38  = simpson38(f, a, b, 99, trueVal);
    auto rRomb = romberg(f, a, b, 6, 1e-10, trueVal);
    auto rG3   = gaussianQuadrature(f, a, b, 3, trueVal);
    auto rG5   = gaussianQuadrature(f, a, b, 5, trueVal);

    std::ostringstream ss;
    ss << "========================================================================================\n";
    ss << "                         NUMERICAL INTEGRATION COMPARISON                               \n";
    ss << "========================================================================================\n";
    if (std::abs(trueVal) > 1e-15) {
        ss << "Reference True Value: " << std::setprecision(10) << trueVal << "\n";
    }
    ss << std::left << std::setw(26) << "Method"
       << std::setw(18) << "Integral Value"
       << std::setw(16) << "Abs Error"
       << std::setw(14) << "Rel Error"
       << std::setw(12) << "Evals" << "\n";
    ss << std::string(86, '-') << "\n";

    auto printRow = [&ss](const IntegrationResult& r) {
        ss << std::left << std::setw(26) << r.method
           << std::setprecision(8) << std::fixed << std::setw(18) << r.value
           << std::scientific << std::setprecision(3) << std::setw(16) << r.absError
           << std::setw(14) << r.relError
           << std::setw(12) << r.functionEvaluations << "\n";
    };

    printRow(rTrap);
    printRow(rS13);
    printRow(rS38);
    printRow(rRomb);
    printRow(rG3);
    printRow(rG5);
    ss << std::string(86, '-') << "\n";

    ss << "QUADRATURE ERROR OBSERVATIONS:\n";
    ss << "- Trapezoidal has global truncation error O(h^2).\n";
    ss << "- Simpson's 1/3 and 3/8 achieve O(h^4) by fitting parabolic/cubic arcs.\n";
    ss << "- Romberg accelerates convergence via Richardson extrapolation, achieving high precision.\n";
    ss << "- Gauss-Legendre quadrature is exact for polynomials up to degree 2n-1 using only n points!\n";

    return ss.str();
}

} // namespace ees::unit2
