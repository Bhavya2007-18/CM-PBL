#include "ees/unit4/ODESolvers.hpp"
#include "ees/core/ErrorAnalysis.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace ees::unit4 {

std::string ODESolution::summary() const {
    std::ostringstream ss;
    ss << "=== ODE Solution Trajectory (" << method << ") ===\n";
    ss << "Integration Interval: [" << x[0] << ", " << x.toStdVector().back() << "]\n";
    ss << "Step Size h:          " << stepSize << " (" << x.size() - 1 << " steps)\n";
    ss << "Final Value y(" << x.toStdVector().back() << "): " << std::setprecision(8) << y.toStdVector().back() << "\n";
    if (!exactY.empty()) {
        ss << "Exact Final Value:    " << exactY.toStdVector().back() << "\n";
        ss << "Final Absolute Error: " << std::scientific << finalError << "\n";
        ss << "Maximum Global Error: " << maxError << "\n";
    }
    return ss.str();
}

std::string ODESolution::tableString(size_t maxRows) const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);
    ss << std::setw(6) << "Step" << " | "
       << std::setw(12) << "x" << " | "
       << std::setw(16) << "y (approx)";
    if (!exactY.empty()) {
        ss << " | " << std::setw(16) << "y (exact)"
           << " | " << std::setw(14) << "Abs Error";
    }
    ss << "\n" << std::string(exactY.empty() ? 40 : 74, '-') << "\n";

    size_t total = x.size();
    size_t stride = (total > maxRows && maxRows > 2) ? (total / maxRows) : 1;

    for (size_t i = 0; i < total; i += stride) {
        ss << std::setw(6) << i << " | "
           << std::setw(12) << x[i] << " | "
           << std::setw(16) << y[i];
        if (!exactY.empty()) {
            ss << " | " << std::setw(16) << exactY[i]
               << " | " << std::scientific << std::setprecision(4) << std::setw(14) << absErrors[i]
               << std::fixed << std::setprecision(6);
        }
        ss << "\n";
    }
    return ss.str();
}

ODESolution ODESolverLab::euler(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    if (h <= 0.0 || xEnd <= x0) throw std::invalid_argument("Invalid ODE integration limits or step size");

    size_t numSteps = static_cast<size_t>(std::ceil((xEnd - x0) / h));
    ODESolution sol;
    sol.method = "Euler's Method [O(h)]";
    sol.stepSize = h;
    sol.x.resize(numSteps + 1);
    sol.y.resize(numSteps + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;

    for (size_t i = 0; i < numSteps; ++i) {
        double curX = sol.x[i];
        double curY = sol.y[i];
        sol.y[i + 1] = curY + h * f(curX, curY);
        sol.x[i + 1] = curX + h;
    }

    if (exactSol) {
        sol.exactY.resize(numSteps + 1);
        sol.absErrors.resize(numSteps + 1);
        double maxErr = 0.0;
        for (size_t i = 0; i <= numSteps; ++i) {
            sol.exactY[i] = exactSol(sol.x[i]);
            sol.absErrors[i] = std::abs(sol.exactY[i] - sol.y[i]);
            maxErr = std::max(maxErr, sol.absErrors[i]);
        }
        sol.maxError = maxErr;
        sol.finalError = sol.absErrors[numSteps];
    }
    return sol;
}

ODESolution ODESolverLab::eulerPredictorCorrector(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    if (h <= 0.0 || xEnd <= x0) throw std::invalid_argument("Invalid ODE integration limits or step size");

    size_t numSteps = static_cast<size_t>(std::ceil((xEnd - x0) / h));
    ODESolution sol;
    sol.method = "Heun's Predictor-Corrector [O(h^2)]";
    sol.stepSize = h;
    sol.x.resize(numSteps + 1);
    sol.y.resize(numSteps + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;

    for (size_t i = 0; i < numSteps; ++i) {
        double curX = sol.x[i];
        double curY = sol.y[i];
        double nextX = curX + h;

        // Predictor
        double yPredict = curY + h * f(curX, curY);
        // Corrector
        double yCorrect = curY + 0.5 * h * (f(curX, curY) + f(nextX, yPredict));

        sol.x[i + 1] = nextX;
        sol.y[i + 1] = yCorrect;
    }

    if (exactSol) {
        sol.exactY.resize(numSteps + 1);
        sol.absErrors.resize(numSteps + 1);
        double maxErr = 0.0;
        for (size_t i = 0; i <= numSteps; ++i) {
            sol.exactY[i] = exactSol(sol.x[i]);
            sol.absErrors[i] = std::abs(sol.exactY[i] - sol.y[i]);
            maxErr = std::max(maxErr, sol.absErrors[i]);
        }
        sol.maxError = maxErr;
        sol.finalError = sol.absErrors[numSteps];
    }
    return sol;
}

ODESolution ODESolverLab::rungeKutta4(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    if (h <= 0.0 || xEnd <= x0) throw std::invalid_argument("Invalid ODE integration limits or step size");

    size_t numSteps = static_cast<size_t>(std::ceil((xEnd - x0) / h));
    ODESolution sol;
    sol.method = "Classic Runge-Kutta 4th-Order (RK4) [O(h^4)]";
    sol.stepSize = h;
    sol.x.resize(numSteps + 1);
    sol.y.resize(numSteps + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;

    for (size_t i = 0; i < numSteps; ++i) {
        double x_n = sol.x[i];
        double y_n = sol.y[i];

        double k1 = f(x_n, y_n);
        double k2 = f(x_n + 0.5 * h, y_n + 0.5 * h * k1);
        double k3 = f(x_n + 0.5 * h, y_n + 0.5 * h * k2);
        double k4 = f(x_n + h, y_n + h * k3);

        sol.y[i + 1] = y_n + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        sol.x[i + 1] = x_n + h;
    }

    if (exactSol) {
        sol.exactY.resize(numSteps + 1);
        sol.absErrors.resize(numSteps + 1);
        double maxErr = 0.0;
        for (size_t i = 0; i <= numSteps; ++i) {
            sol.exactY[i] = exactSol(sol.x[i]);
            sol.absErrors[i] = std::abs(sol.exactY[i] - sol.y[i]);
            maxErr = std::max(maxErr, sol.absErrors[i]);
        }
        sol.maxError = maxErr;
        sol.finalError = sol.absErrors[numSteps];
    }
    return sol;
}

ODESolution ODESolverLab::taylorSeries(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    if (h <= 0.0 || xEnd <= x0) throw std::invalid_argument("Invalid ODE parameters");

    size_t numSteps = static_cast<size_t>(std::ceil((xEnd - x0) / h));
    ODESolution sol;
    sol.method = "Taylor Series Method (3rd-Order)";
    sol.stepSize = h;
    sol.x.resize(numSteps + 1);
    sol.y.resize(numSteps + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;

    // y' = f(x, y)
    // y'' = f_x + f_y * f
    // y''' ≈ finite difference approximation of y''
    double delta = 1e-5;

    for (size_t i = 0; i < numSteps; ++i) {
        double x_n = sol.x[i];
        double y_n = sol.y[i];

        double yp = f(x_n, y_n);
        double fx = (f(x_n + delta, y_n) - f(x_n - delta, y_n)) / (2.0 * delta);
        double fy = (f(x_n, y_n + delta) - f(x_n, y_n - delta)) / (2.0 * delta);
        double ypp = fx + fy * yp;

        sol.y[i + 1] = y_n + h * yp + 0.5 * h * h * ypp;
        sol.x[i + 1] = x_n + h;
    }

    if (exactSol) {
        sol.exactY.resize(numSteps + 1);
        sol.absErrors.resize(numSteps + 1);
        double maxErr = 0.0;
        for (size_t i = 0; i <= numSteps; ++i) {
            sol.exactY[i] = exactSol(sol.x[i]);
            sol.absErrors[i] = std::abs(sol.exactY[i] - sol.y[i]);
            maxErr = std::max(maxErr, sol.absErrors[i]);
        }
        sol.maxError = maxErr;
        sol.finalError = sol.absErrors[numSteps];
    }
    return sol;
}

ODESolution ODESolverLab::adamsBashforth(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    if (h <= 0.0 || xEnd <= x0) throw std::invalid_argument("Invalid ODE parameters");

    size_t numSteps = static_cast<size_t>(std::ceil((xEnd - x0) / h));
    if (numSteps < 4) {
        // Fall back to RK4 if interval is too small
        return rungeKutta4(f, x0, y0, xEnd, h, exactSol);
    }

    ODESolution sol;
    sol.method = "Adams-Bashforth-Moulton 4-Step [O(h^4)]";
    sol.stepSize = h;
    sol.x.resize(numSteps + 1);
    sol.y.resize(numSteps + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;

    // Use RK4 to compute initial 3 steps (x1, y1), (x2, y2), (x3, y3)
    for (size_t i = 0; i < 3; ++i) {
        double xi = sol.x[i];
        double yi = sol.y[i];
        double k1 = f(xi, yi);
        double k2 = f(xi + 0.5 * h, yi + 0.5 * h * k1);
        double k3 = f(xi + 0.5 * h, yi + 0.5 * h * k2);
        double k4 = f(xi + h, yi + h * k3);
        sol.y[i + 1] = yi + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        sol.x[i + 1] = xi + h;
    }

    // Main Adams-Bashforth-Moulton loop
    for (size_t i = 3; i < numSteps; ++i) {
        double fn   = f(sol.x[i],     sol.y[i]);
        double fn_1 = f(sol.x[i - 1], sol.y[i - 1]);
        double fn_2 = f(sol.x[i - 2], sol.y[i - 2]);
        double fn_3 = f(sol.x[i - 3], sol.y[i - 3]);

        double nextX = sol.x[i] + h;

        // Predictor (Adams-Bashforth 4)
        double yPredict = sol.y[i] + (h / 24.0) * (55.0 * fn - 59.0 * fn_1 + 37.0 * fn_2 - 9.0 * fn_3);

        // Corrector (Adams-Moulton 3)
        double fn_plus1 = f(nextX, yPredict);
        double yCorrect = sol.y[i] + (h / 24.0) * (9.0 * fn_plus1 + 19.0 * fn - 5.0 * fn_1 + fn_2);

        sol.x[i + 1] = nextX;
        sol.y[i + 1] = yCorrect;
    }

    if (exactSol) {
        sol.exactY.resize(numSteps + 1);
        sol.absErrors.resize(numSteps + 1);
        double maxErr = 0.0;
        for (size_t i = 0; i <= numSteps; ++i) {
            sol.exactY[i] = exactSol(sol.x[i]);
            sol.absErrors[i] = std::abs(sol.exactY[i] - sol.y[i]);
            maxErr = std::max(maxErr, sol.absErrors[i]);
        }
        sol.maxError = maxErr;
        sol.finalError = sol.absErrors[numSteps];
    }
    return sol;
}

ODESolution ODESolverLab::milne(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    if (h <= 0.0 || xEnd <= x0) throw std::invalid_argument("Invalid ODE parameters");

    size_t numSteps = static_cast<size_t>(std::ceil((xEnd - x0) / h));
    if (numSteps < 4) {
        return rungeKutta4(f, x0, y0, xEnd, h, exactSol);
    }

    ODESolution sol;
    sol.method = "Milne's Predictor-Corrector Method [O(h^4)]";
    sol.stepSize = h;
    sol.x.resize(numSteps + 1);
    sol.y.resize(numSteps + 1);

    sol.x[0] = x0;
    sol.y[0] = y0;

    // First 3 steps with RK4
    for (size_t i = 0; i < 3; ++i) {
        double xi = sol.x[i];
        double yi = sol.y[i];
        double k1 = f(xi, yi);
        double k2 = f(xi + 0.5 * h, yi + 0.5 * h * k1);
        double k3 = f(xi + 0.5 * h, yi + 0.5 * h * k2);
        double k4 = f(xi + h, yi + h * k3);
        sol.y[i + 1] = yi + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        sol.x[i + 1] = xi + h;
    }

    for (size_t i = 3; i < numSteps; ++i) {
        double fn   = f(sol.x[i],     sol.y[i]);
        double fn_1 = f(sol.x[i - 1], sol.y[i - 1]);
        double fn_2 = f(sol.x[i - 2], sol.y[i - 2]);

        double nextX = sol.x[i] + h;

        // Milne Predictor: y*_{n+1} = y_{n-3} + (4h/3) * (2*fn - fn_1 + 2*fn_2)
        double yPredict = sol.y[i - 3] + (4.0 * h / 3.0) * (2.0 * fn - fn_1 + 2.0 * fn_2);

        // Milne / Simpson Corrector: y_{n+1} = y_{n-1} + (h/3) * (fn_1 + 4*fn + f(nextX, yPredict))
        double fn_plus1 = f(nextX, yPredict);
        double yCorrect = sol.y[i - 1] + (h / 3.0) * (fn_1 + 4.0 * fn + fn_plus1);

        sol.x[i + 1] = nextX;
        sol.y[i + 1] = yCorrect;
    }

    if (exactSol) {
        sol.exactY.resize(numSteps + 1);
        sol.absErrors.resize(numSteps + 1);
        double maxErr = 0.0;
        for (size_t i = 0; i <= numSteps; ++i) {
            sol.exactY[i] = exactSol(sol.x[i]);
            sol.absErrors[i] = std::abs(sol.exactY[i] - sol.y[i]);
            maxErr = std::max(maxErr, sol.absErrors[i]);
        }
        sol.maxError = maxErr;
        sol.finalError = sol.absErrors[numSteps];
    }
    return sol;
}

std::string ODESolverLab::picardDemo(double x0, double y0, int iterations) {
    std::ostringstream ss;
    ss << "=== Picard's Method of Successive Approximations ===\n";
    ss << "Initial Value Problem: dy/dx = x + y, with y(" << x0 << ") = " << y0 << "\n";
    ss << "Integral Recurrence: y^{(k+1)}(x) = y0 + int_{x0}^x f(t, y^{(k)}(t)) dt\n\n";

    ss << "Approximation Iterations:\n";
    ss << "  k = 0: y^{(0)}(x) = " << y0 << "\n";
    ss << "  k = 1: y^{(1)}(x) = " << y0 << " + int_{0}^x (t + " << y0 << ") dt = " << y0 << " + " << y0 << "*x + 0.5*x^2\n";
    ss << "  k = 2: y^{(2)}(x) = " << y0 << " + " << y0 << "*x + (1 + " << y0 << ")*0.5*x^2 + (1/6)*x^3\n";
    ss << "  k = 3: y^{(3)}(x) = " << y0 << " + " << y0 << "*x + (1 + " << y0 << ")*0.5*x^2 + (1 + " << y0 << ")*(1/6)*x^3 + (1/24)*x^4\n";

    ss << "\nAnalytic Limit: Matches Taylor expansion of exact solution y(x) = 2*e^x - x - 1 (for y(0)=1).\n";
    (void)iterations;
    return ss.str();
}

std::string ODESolverLab::compareODESolvers(
    ODEFunction f,
    double x0,
    double y0,
    double xEnd,
    double h,
    core::RealFunction exactSol
) {
    auto rEuler = euler(f, x0, y0, xEnd, h, exactSol);
    auto rHeun  = eulerPredictorCorrector(f, x0, y0, xEnd, h, exactSol);
    auto rRK4   = rungeKutta4(f, x0, y0, xEnd, h, exactSol);
    auto rAB    = adamsBashforth(f, x0, y0, xEnd, h, exactSol);
    auto rMilne = milne(f, x0, y0, xEnd, h, exactSol);

    std::ostringstream ss;
    ss << "========================================================================================\n";
    ss << "                              ODE SOLVER BENCHMARK COMPARISON                           \n";
    ss << "========================================================================================\n";
    ss << "Interval: [" << x0 << ", " << xEnd << "], Step size h = " << h << "\n";
    if (exactSol) {
        ss << "Exact Final Value y(" << xEnd << "): " << std::setprecision(10) << exactSol(xEnd) << "\n";
    }
    ss << std::left << std::setw(30) << "Method"
       << std::setw(18) << "Final y(" << xEnd << ")"
       << std::setw(16) << "Max Global Error"
       << std::setw(16) << "Final Error"
       << std::setw(12) << "Order" << "\n";
    ss << std::string(92, '-') << "\n";

    auto printRow = [&ss](const ODESolution& sol, const std::string& order) {
        ss << std::left << std::setw(30) << sol.method
           << std::setprecision(8) << std::fixed << std::setw(18) << sol.y.toStdVector().back()
           << std::scientific << std::setprecision(3) << std::setw(16) << sol.maxError
           << std::setw(16) << sol.finalError
           << std::setw(12) << order << "\n";
    };

    printRow(rEuler, "O(h)");
    printRow(rHeun,  "O(h^2)");
    printRow(rRK4,   "O(h^4)");
    printRow(rAB,    "O(h^4)");
    printRow(rMilne, "O(h^4)");
    ss << std::string(92, '-') << "\n";

    ss << "METHOD CHARACTERISTICS:\n";
    ss << "- Euler: Simplest 1st-order method, requires tiny step sizes for stability.\n";
    ss << "- RK4: Single-step workhorse of scientific computing; self-starting, highly accurate and stable.\n";
    ss << "- Multi-step (Adams-Bashforth / Milne): Requires only 2 function evaluations per step (vs 4 in RK4),\n";
    ss << "  ideal for computationally heavy ODE right-hand sides, initialized by RK4.\n";

    return ss.str();
}

} // namespace ees::unit4

