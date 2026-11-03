#include "ees/unit3/Splines.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace ees::unit3 {

std::string SplineResult::summary() const {
    std::ostringstream ss;
    ss << "=== Spline Approximation (" << splineType << ") ===\n";
    ss << "Evaluated at x = " << std::fixed << std::setprecision(6) << evaluatedX
       << " => y = " << interpolatedY << "\n";

    if (!segments.empty()) {
        ss << "\nCubic Spline Piecewise Polynomials:\n";
        for (size_t i = 0; i < segments.size(); ++i) {
            const auto& s = segments[i];
            ss << "Interval " << (i + 1) << " [x" << i << "]: S_" << i << "(x) = "
               << s.a << " + (" << s.b << ")*(x - " << s.x_i << ") + ("
               << s.c << ")*(x - " << s.x_i << ")^2 + ("
               << s.d << ")*(x - " << s.x_i << ")^3\n";
        }
    }
    return ss.str();
}

SplineResult SplineLab::linearSpline(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    const size_t n = x.size();
    if (n < 2 || n != y.size()) throw std::invalid_argument("Linear spline requires >= 2 points");

    SplineResult res;
    res.splineType = "First-Degree (Linear) Spline";
    res.evaluatedX = xVal;

    // Boundary clamp or locate segment
    size_t idx = 0;
    if (xVal <= x[0]) {
        idx = 0;
    } else if (xVal >= x[n - 1]) {
        idx = n - 2;
    } else {
        while (idx < n - 2 && xVal > x[idx + 1]) {
            idx++;
        }
    }

    double slope = (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]);
    res.interpolatedY = y[idx] + slope * (xVal - x[idx]);
    return res;
}

SplineResult SplineLab::quadraticSpline(
    const core::Vector& x,
    const core::Vector& y,
    double xVal,
    double initialM0
) {
    const size_t n = x.size();
    if (n < 3 || n != y.size()) throw std::invalid_argument("Quadratic spline requires >= 3 points");

    // S_i(x) = a_i * (x - x_i)^2 + b_i * (x - x_i) + c_i
    // With c_i = y_i, b_i = m_i (derivative at x_i), and a_i = (m_{i+1} - m_i) / (2 * h_i)
    // Continuity of y gives: y_{i+1} = a_i * h_i^2 + b_i * h_i + c_i
    // => m_{i+1} = 2 * (y_{i+1} - y_i)/h_i - m_i
    std::vector<double> m(n);
    m[0] = initialM0;

    for (size_t i = 0; i < n - 1; ++i) {
        double hi = x[i + 1] - x[i];
        m[i + 1] = 2.0 * (y[i + 1] - y[i]) / hi - m[i];
    }

    SplineResult res;
    res.splineType = "Second-Degree (Quadratic) Spline";
    res.evaluatedX = xVal;

    size_t idx = 0;
    if (xVal <= x[0]) idx = 0;
    else if (xVal >= x[n - 1]) idx = n - 2;
    else {
        while (idx < n - 2 && xVal > x[idx + 1]) idx++;
    }

    double hi = x[idx + 1] - x[idx];
    double ai = (m[idx + 1] - m[idx]) / (2.0 * hi);
    double bi = m[idx];
    double ci = y[idx];
    double dx = xVal - x[idx];

    res.interpolatedY = ai * dx * dx + bi * dx + ci;
    return res;
}

SplineResult SplineLab::naturalCubicSpline(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    const size_t n = x.size();
    if (n < 3 || n != y.size()) throw std::invalid_argument("Natural cubic spline requires >= 3 points");

    // Number of intervals = n - 1
    const size_t numInt = n - 1;
    std::vector<double> h(numInt);
    for (size_t i = 0; i < numInt; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    // Set up tridiagonal system for second derivatives M_i = S''(x_i):
    // h_{i-1}*M_{i-1} + 2*(h_{i-1} + h_i)*M_i + h_i*M_{i+1} = 6 * [(y_{i+1}-y_i)/h_i - (y_i - y_{i-1})/h_{i-1}]
    // Natural boundary condition: M_0 = 0, M_{n-1} = 0
    std::vector<double> alpha(n, 0.0);
    for (size_t i = 1; i < numInt; ++i) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    // Solve tridiagonal system via Thomas algorithm
    std::vector<double> l(n, 1.0);
    std::vector<double> mu(n, 0.0);
    std::vector<double> z(n, 0.0);

    for (size_t i = 1; i < numInt; ++i) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    std::vector<double> c(n, 0.0); // c_i = M_i / 2
    std::vector<double> b(numInt, 0.0);
    std::vector<double> d(numInt, 0.0);

    // Back substitution
    for (int j = static_cast<int>(numInt) - 1; j >= 1; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
    }

    SplineResult res;
    res.splineType = "Natural Cubic Spline (C2 Continuous)";
    res.evaluatedX = xVal;

    for (size_t i = 0; i < numInt; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);

        CubicSplineSegment seg;
        seg.x_i = x[i];
        seg.a = y[i];
        seg.b = b[i];
        seg.c = c[i];
        seg.d = d[i];
        res.segments.push_back(seg);
    }

    // Evaluate at target xVal
    size_t idx = 0;
    if (xVal <= x[0]) idx = 0;
    else if (xVal >= x[n - 1]) idx = numInt - 1;
    else {
        while (idx < numInt - 1 && xVal > x[idx + 1]) idx++;
    }

    double dx = xVal - x[idx];
    res.interpolatedY = y[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
    return res;
}

core::Vector SplineLab::evaluateNaturalCubicSpline(
    const core::Vector& xKnots,
    const core::Vector& yKnots,
    const core::Vector& evalPoints
) {
    core::Vector result(evalPoints.size());
    for (size_t i = 0; i < evalPoints.size(); ++i) {
        result[i] = naturalCubicSpline(xKnots, yKnots, evalPoints[i]).interpolatedY;
    }
    return result;
}

core::Vector SplineLab::evaluateBSpline(
    const core::Vector& controlPointsX,
    const core::Vector& controlPointsY,
    int numSamples,
    int degree
) {
    // Generates B-Spline curve samples using Cox-de Boor recursive formulation
    (void)degree;
    size_t n = controlPointsX.size();
    if (n < 4) throw std::invalid_argument("Cubic B-spline requires at least 4 control points");

    core::Vector curveY(numSamples);
    double step = static_cast<double>(n - 3) / static_cast<double>(numSamples - 1);

    for (int s = 0; s < numSamples; ++s) {
        double u = s * step;
        int i = static_cast<int>(std::floor(u));
        if (i >= static_cast<int>(n) - 3) i = static_cast<int>(n) - 4;
        double t = u - i; // Local parameter in [0, 1]

        // Uniform cubic B-spline blending basis matrix:
        // [1 t t^2 t^3] * (1/6) * [ 1  4  1  0 ]
        //                         [-3  0  3  0 ]
        //                         [ 3 -6  3  0 ]
        //                         [-1  3 -3  1 ]
        double B0 = (1.0 - t) * (1.0 - t) * (1.0 - t) / 6.0;
        double B1 = (3.0 * t * t * t - 6.0 * t * t + 4.0) / 6.0;
        double B2 = (-3.0 * t * t * t + 3.0 * t * t + 3.0 * t + 1.0) / 6.0;
        double B3 = (t * t * t) / 6.0;

        curveY[s] = B0 * controlPointsY[i] + B1 * controlPointsY[i + 1] +
                    B2 * controlPointsY[i + 2] + B3 * controlPointsY[i + 3];
    }
    return curveY;
}

std::string SplineLab::compareSplines(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    auto lin = linearSpline(x, y, xVal);
    auto quad = quadraticSpline(x, y, xVal);
    auto cubic = naturalCubicSpline(x, y, xVal);

    std::ostringstream ss;
    ss << "========================================================================================\n";
    ss << "                          SPLINE APPROXIMATION COMPARISON                               \n";
    ss << "========================================================================================\n";
    ss << "Evaluation Point: x = " << xVal << "\n";
    ss << std::left << std::setw(28) << "Spline Degree"
       << std::setw(20) << "Interpolated y"
       << std::setw(35) << "Continuity Class" << "\n";
    ss << std::string(83, '-') << "\n";

    ss << std::left << std::setw(28) << "First-Degree (Linear)"
       << std::setprecision(8) << std::fixed << std::setw(20) << lin.interpolatedY
       << std::setw(35) << "C^0 (Continuous, piecewise linear)" << "\n";

    ss << std::left << std::setw(28) << "Second-Degree (Quadratic)"
       << std::setprecision(8) << std::fixed << std::setw(20) << quad.interpolatedY
       << std::setw(35) << "C^1 (Continuous first derivative)" << "\n";

    ss << std::left << std::setw(28) << "Third-Degree (Natural Cubic)"
       << std::setprecision(8) << std::fixed << std::setw(20) << cubic.interpolatedY
       << std::setw(35) << "C^2 (Continuous curvature, minimal energy)" << "\n";
    ss << std::string(83, '-') << "\n";

    ss << "THEORETICAL INSIGHT (Interpolation vs Approximation):\n";
    ss << "- Interpolation splines (Natural Cubic) pass EXACTLY through all knot points.\n";
    ss << "- B-Splines provide smooth approximation within the convex hull of control points,\n";
    ss << "  preventing Runge's oscillation phenomenon common in high-degree global polynomials.\n";

    return ss.str();
}

} // namespace ees::unit3

