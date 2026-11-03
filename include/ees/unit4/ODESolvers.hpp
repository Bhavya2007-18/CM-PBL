#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/NumericalResult.hpp"
#include "ees/core/ExpressionParser.hpp"
#include <vector>
#include <string>

namespace ees::unit4 {

// dy/dx = f(x, y)
using ODEFunction = std::function<double(double, double)>;

struct ODESolution {
    core::Vector x;
    core::Vector y;
    core::Vector exactY;
    core::Vector absErrors;
    double maxError{0.0};
    double finalError{0.0};
    double stepSize{0.0};
    std::string method;
    std::string summary() const;
    std::string tableString(size_t maxRows = 20) const;
};

class ODESolverLab {
public:
    // Euler's Method: y_{n+1} = y_n + h * f(x_n, y_n), Global error O(h)
    static ODESolution euler(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );

    // Heun's Method (Euler Predictor-Corrector / Improved Euler): Global error O(h^2)
    // Predictor: y*_{n+1} = y_n + h * f(x_n, y_n)
    // Corrector: y_{n+1} = y_n + (h/2) * [f(x_n, y_n) + f(x_{n+1}, y*_{n+1})]
    static ODESolution eulerPredictorCorrector(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );

    // Classic 4th-Order Runge-Kutta Method (RK4): Global error O(h^4)
    // k1 = f(x, y)
    // k2 = f(x + h/2, y + h*k1/2)
    // k3 = f(x + h/2, y + h*k2/2)
    // k4 = f(x + h, y + h*k3)
    // y_{n+1} = y_n + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    static ODESolution rungeKutta4(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );

    // Taylor Series Method (up to 3 terms using numerical derivatives): Global error O(h^3)
    // y_{n+1} = y_n + h*y' + (h^2/2)*y'' + (h^3/6)*y'''
    static ODESolution taylorSeries(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );

    // Adams-Bashforth-Moulton 4-step Predictor-Corrector: Global error O(h^4)
    // Starts with RK4 for first 3 steps, then:
    // Predictor (Adams-Bashforth 4): y*_{n+1} = y_n + (h/24) * [55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3}]
    // Corrector (Adams-Moulton 3):   y_{n+1} = y_n + (h/24) * [9 f_{n+1} + 19 f_n - 5 f_{n-1} + f_{n-2}]
    static ODESolution adamsBashforth(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );

    // Milne's Predictor-Corrector Method: Global error O(h^4)
    // Starts with RK4 for first 3 steps, then:
    // Predictor: y*_{n+1} = y_{n-3} + (4h/3) * [2 f_n - f_{n-1} + 2 f_{n-2}]
    // Corrector: y_{n+1} = y_{n-1} + (h/3) * [f_{n-1} + 4 f_n + f(x_{n+1}, y*_{n+1})]
    static ODESolution milne(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );

    // Picard's Successive Approximations demonstration:
    // y^{(k+1)}(x) = y0 + int_{x0}^x f(t, y^{(k)}(t)) dt
    static std::string picardDemo(
        double x0,
        double y0,
        int iterations = 4
    );

    // Compare all ODE solvers on the same initial value problem
    static std::string compareODESolvers(
        ODEFunction f,
        double x0,
        double y0,
        double xEnd,
        double h,
        core::RealFunction exactSol = nullptr
    );
};

} // namespace ees::unit4

