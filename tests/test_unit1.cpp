#include "ees/unit1/Foundations.hpp"
#include "ees/unit1/RootFinding.hpp"
#include "ees/unit1/Optimization.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace ees::unit1;

void testFoundations() {
    std::cout << "[TEST] Running Unit I Foundations & Taylor Series tests..." << std::endl;

    // f(x) = exp(x), x0 = 0, evaluate at x = 1 (approx e ≈ 2.7182818)
    auto fExp = [](double x){ return std::exp(x); };
    auto taylorRes = FoundationsLab::computeTaylorSeries(fExp, 0.0, 1.0, 5, std::exp(1.0));
    assert(taylorRes.terms.size() == 6);
    // Taylor for exp(1) with 5 orders: 1 + 1 + 1/2 + 1/6 + 1/24 + 1/120 = 2.716666...
    assert(std::abs(taylorRes.approxValue - 2.716666666) < 1e-4);

    // Rolle's Theorem on f(x) = x*(x - 2) on [0, 2]: f(0)=0, f(2)=0. f'(c) = 2c - 2 = 0 => c = 1.0
    auto fParabola = [](double x){ return x * (x - 2.0); };
    auto rolleRes = FoundationsLab::verifyRollesTheorem(fParabola, 0.0, 2.0);
    assert(rolleRes.theoremHolds);
    assert(!rolleRes.cValues.empty());
    assert(std::abs(rolleRes.cValues[0] - 1.0) < 1e-3);

    // Mean Value Theorem on f(x) = x^3 on [0, 3]: (27 - 0)/3 = 9. f'(c) = 3c^2 = 9 => c = sqrt(3) ≈ 1.73205
    auto fCubic = [](double x){ return x * x * x; };
    auto mvtRes = FoundationsLab::verifyMeanValueTheorem(fCubic, 0.0, 3.0);
    assert(mvtRes.theoremHolds);
    assert(!mvtRes.cValues.empty());
    assert(std::abs(mvtRes.cValues[0] - std::sqrt(3.0)) < 1e-3);

    std::cout << "  -> Unit I Foundations tests PASSED." << std::endl;
}

void testRootFinding() {
    std::cout << "[TEST] Running Unit I Root Finding tests..." << std::endl;

    // Problem: x^3 - 4x - 9 = 0. Known root in [2, 3] is r ≈ 2.70652796
    auto f = [](double x){ return x * x * x - 4.0 * x - 9.0; };
    const double expectedRoot = 2.70652796;

    // 1. Bisection
    auto resBis = RootFindingLab::bisection(f, 2.0, 3.0, 1e-6);
    assert(resBis.success);
    assert(std::abs(resBis.solution - expectedRoot) < 1e-5);
    assert(resBis.convergenceRate == ees::core::ConvergenceRate::Linear);

    // 2. Newton-Raphson
    auto resNew = RootFindingLab::newtonRaphson(f, 2.5, 1e-7);
    assert(resNew.success);
    assert(std::abs(resNew.solution - expectedRoot) < 1e-6);
    assert(resNew.iterations < resBis.iterations); // Quadratic is faster
    assert(resNew.convergenceRate == ees::core::ConvergenceRate::Quadratic);

    // 3. Secant
    auto resSec = RootFindingLab::secant(f, 2.0, 3.0, 1e-7);
    assert(resSec.success);
    assert(std::abs(resSec.solution - expectedRoot) < 1e-6);
    assert(resSec.convergenceRate == ees::core::ConvergenceRate::Superlinear);

    std::cout << "  -> Unit I Root Finding tests PASSED." << std::endl;
}

void testOptimization() {
    std::cout << "[TEST] Running Unit I Optimization tests..." << std::endl;

    // 1D Minima: f(x) = (x - 3)^2 + 4. Minimum at x = 3, min value = 4.
    auto f1 = [](double x){ return (x - 3.0) * (x - 3.0) + 4.0; };

    auto fibRes = OptimizationLab::fibonacciSearch(f1, 0.0, 5.0, 1e-4);
    assert(fibRes.success);
    assert(std::abs(fibRes.solution - 3.0) < 1e-3);

    auto goldRes = OptimizationLab::goldenSectionSearch(f1, 0.0, 5.0, 1e-4);
    assert(goldRes.success);
    assert(std::abs(goldRes.solution - 3.0) < 1e-3);

    auto newtRes = OptimizationLab::newtonsMethod1D(f1, 1.0, 1e-6);
    assert(newtRes.success);
    assert(std::abs(newtRes.solution - 3.0) < 1e-5);

    // Multivariate: Rosenbrock function or sphere f(x, y) = (x-2)^2 + (y-3)^2
    // Minimum at (2, 3), min value = 0.
    auto fMulti = [](const ees::core::Vector& v){
        return (v[0] - 2.0) * (v[0] - 2.0) + (v[1] - 3.0) * (v[1] - 3.0);
    };

    ees::core::Vector start{0.0, 0.0};

    // Steepest Descent
    auto sdRes = OptimizationLab::steepestDescent(fMulti, start, 0.2, 1e-4);
    assert(sdRes.success);
    assert(std::abs(sdRes.solutionVector[0] - 2.0) < 1e-3);
    assert(std::abs(sdRes.solutionVector[1] - 3.0) < 1e-3);

    // Nelder-Mead
    auto nmRes = OptimizationLab::nelderMead(fMulti, start, 0.5, 1e-4);
    assert(nmRes.success);
    assert(std::abs(nmRes.solutionVector[0] - 2.0) < 1e-3);
    assert(std::abs(nmRes.solutionVector[1] - 3.0) < 1e-3);

    std::cout << "  -> Unit I Optimization tests PASSED." << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "    EES UNIT I AUTOMATED TEST SUITE     " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testFoundations();
        testRootFinding();
        testOptimization();

        std::cout << "\n========================================" << std::endl;
        std::cout << "    ALL UNIT I TESTS PASSED (3/3)       " << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "UNIT I TEST FAILURE: " << ex.what() << std::endl;
        return 1;
    }
}

