#include "ees/unit4/ODESolvers.hpp"
#include "ees/unit4/PDESolvers.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace ees::unit4;

void testODESolvers() {
    std::cout << "[TEST] Running Unit IV ODE Solvers tests..." << std::endl;

    // IVP: dy/dx = y, y(0) = 1. Exact solution: y(x) = exp(x).
    auto f = [](double x, double y) {
        (void)x;
        return y;
    };
    auto exact = [](double x) { return std::exp(x); };

    double x0 = 0.0, y0 = 1.0, xEnd = 1.0, h = 0.1;
    double expectedFinal = std::exp(1.0); // 2.718281828...

    // 1. Euler
    auto rEuler = ODESolverLab::euler(f, x0, y0, xEnd, h, exact);
    assert(rEuler.y.size() == 11);
    assert(std::abs(rEuler.y.toStdVector().back() - expectedFinal) < 0.2); // Euler error O(h)

    // 2. Heun (Predictor-Corrector)
    auto rHeun = ODESolverLab::eulerPredictorCorrector(f, x0, y0, xEnd, h, exact);
    assert(std::abs(rHeun.y.toStdVector().back() - expectedFinal) < 0.01); // Heun error O(h^2)

    // 3. RK4
    auto rRK4 = ODESolverLab::rungeKutta4(f, x0, y0, xEnd, h, exact);
    assert(std::abs(rRK4.y.toStdVector().back() - expectedFinal) < 1e-5); // RK4 error O(h^4)

    // 4. Adams-Bashforth-Moulton
    auto rAB = ODESolverLab::adamsBashforth(f, x0, y0, xEnd, h, exact);
    assert(std::abs(rAB.y.toStdVector().back() - expectedFinal) < 1e-4);

    // 5. Milne's Predictor-Corrector
    auto rMilne = ODESolverLab::milne(f, x0, y0, xEnd, h, exact);
    assert(std::abs(rMilne.y.toStdVector().back() - expectedFinal) < 1e-4);

    // Verify error hierarchy: Error(Euler) > Error(Heun) > Error(RK4)
    assert(rEuler.finalError > rHeun.finalError);
    assert(rHeun.finalError > rRK4.finalError);

    std::cout << "  -> Unit IV ODE Solvers tests PASSED." << std::endl;
}

void testPDESolvers() {
    std::cout << "[TEST] Running Unit IV PDE Solvers tests..." << std::endl;

    // 1. Heat Equation 1D: u(x, 0) = sin(pi * x), ends at 0
    // Analytic solution: u(x, t) = exp(-alpha * pi^2 * t) * sin(pi * x)
    auto u0 = [](double x){ return std::sin(3.141592653589793 * x); };
    double alpha = 0.1;
    auto heatSol = PDESolverLab::solveHeatEquation1D(alpha, 1.0, 0.5, 21, 101, u0, 0.0, 0.0);
    assert(heatSol.isStable);
    assert(heatSol.grid.rows() == 101);
    assert(heatSol.grid.cols() == 21);
    // Center point x = 0.5 starts at 1.0 and decays over time
    double initialCenter = heatSol.grid(0, 10);
    double finalCenter = heatSol.grid(100, 10);
    assert(std::abs(initialCenter - 1.0) < 1e-4);
    assert(finalCenter < initialCenter); // Temperature must diffuse/decay!
    assert(finalCenter > 0.0);

    // 2. Laplace 2D: Top = 100, other edges = 0
    auto laplaceSol = PDESolverLab::solveLaplaceEquation2D(1.0, 1.0, 11, 11, 100.0, 0.0, 0.0, 0.0, 1.3, 1e-4, 500);
    assert(laplaceSol.grid.rows() == 11);
    assert(laplaceSol.grid.cols() == 11);
    // Boundary checks
    assert(std::abs(laplaceSol.grid(10, 5) - 100.0) < 1e-4); // Top
    assert(std::abs(laplaceSol.grid(0, 5) - 0.0) < 1e-4);    // Bottom
    // Interior temperature must be strictly between 0 and 100 (Maximum Principle)
    double centerTemp = laplaceSol.grid(5, 5);
    assert(centerTemp > 0.0 && centerTemp < 100.0);

    // 3. Wave Equation 1D: u(x, 0) = sin(pi * x), du/dt = 0
    // Analytic solution: u(x, t) = sin(pi * x) * cos(pi * c * t)
    // For c = 1, at period T = 2, wave returns to initial state!
    auto waveSol = PDESolverLab::solveWaveEquation1D(1.0, 1.0, 1.0, 21, 101, u0, nullptr, 0.0, 0.0);
    assert(waveSol.isStable);
    assert(waveSol.grid.rows() == 101);

    std::cout << "  -> Unit IV PDE Solvers tests PASSED." << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "    EES UNIT IV AUTOMATED TEST SUITE    " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testODESolvers();
        testPDESolvers();

        std::cout << "\n========================================" << std::endl;
        std::cout << "    ALL UNIT IV TESTS PASSED (2/2)      " << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "UNIT IV TEST FAILURE: " << ex.what() << std::endl;
        return 1;
    }
}

