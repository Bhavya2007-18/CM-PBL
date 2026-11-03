#include "ees/unit2/Interpolation.hpp"
#include "ees/unit2/Integration.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace ees::unit2;

void testInterpolation() {
    std::cout << "[TEST] Running Unit II Interpolation tests..." << std::endl;

    // Test dataset from known polynomial y = x^3 - 2x + 1
    // x = [1, 2, 3, 4]
    // y = [0, 5, 22, 57]
    ees::core::Vector x{1.0, 2.0, 3.0, 4.0};
    ees::core::Vector y{0.0, 5.0, 22.0, 57.0};

    // Forward difference table check
    auto fTable = InterpolationLab::computeForwardDifferenceTable(x, y);
    assert(fTable.isEquispaced);
    assert(std::abs(fTable.h - 1.0) < 1e-12);
    // Delta y0 = 5 - 0 = 5
    assert(std::abs(fTable.diffTable(0, 1) - 5.0) < 1e-12);

    // Evaluate at x = 2.5: true y = 2.5^3 - 2*(2.5) + 1 = 15.625 - 5 + 1 = 11.625
    double targetX = 2.5;
    double expectedY = 11.625;

    // 1. Gregory-Newton Forward
    auto gnF = InterpolationLab::gregoryNewtonForward(x, y, targetX);
    assert(std::abs(gnF.interpolatedY - expectedY) < 1e-6);

    // 2. Gregory-Newton Backward
    auto gnB = InterpolationLab::gregoryNewtonBackward(x, y, targetX);
    assert(std::abs(gnB.interpolatedY - expectedY) < 1e-6);

    // 3. Lagrange
    auto lag = InterpolationLab::lagrange(x, y, targetX);
    assert(std::abs(lag.interpolatedY - expectedY) < 1e-6);

    // 4. Newton Divided Difference
    auto ndd = InterpolationLab::newtonDividedDifference(x, y, targetX);
    assert(std::abs(ndd.interpolatedY - expectedY) < 1e-6);

    // Non-equispaced points test with Newton Divided Difference
    ees::core::Vector xUnequal{1.0, 1.5, 3.0, 5.0};
    ees::core::Vector yUnequal{
        1.0 * 1.0,
        1.5 * 1.5,
        3.0 * 3.0,
        5.0 * 5.0
    }; // y = x^2
    auto nddUnequal = InterpolationLab::newtonDividedDifference(xUnequal, yUnequal, 2.0);
    assert(std::abs(nddUnequal.interpolatedY - 4.0) < 1e-6);

    std::cout << "  -> Unit II Interpolation tests PASSED." << std::endl;
}

void testIntegration() {
    std::cout << "[TEST] Running Unit II Numerical Integration tests..." << std::endl;

    // Integral: int_0^1 (4 / (1 + x^2)) dx = 4 * atan(1) = pi ≈ 3.141592653589793
    auto fPi = [](double x){ return 4.0 / (1.0 + x * x); };
    const double truePi = 3.14159265358979323846;

    // 1. Trapezoidal
    auto resTrap = IntegrationLab::trapezoidal(fPi, 0.0, 1.0, 100, truePi);
    assert(std::abs(resTrap.value - truePi) < 1e-4);

    // 2. Simpson's 1/3
    auto resS13 = IntegrationLab::simpson13(fPi, 0.0, 1.0, 100, truePi);
    assert(std::abs(resS13.value - truePi) < 1e-7);

    // 3. Simpson's 3/8
    auto resS38 = IntegrationLab::simpson38(fPi, 0.0, 1.0, 99, truePi);
    assert(std::abs(resS38.value - truePi) < 1e-7);

    // 4. Romberg
    auto resRomb = IntegrationLab::romberg(fPi, 0.0, 1.0, 6, 1e-10, truePi);
    assert(std::abs(resRomb.value - truePi) < 1e-9);

    // 5. Gaussian Quadrature (3-point and 5-point)
    auto resG3 = IntegrationLab::gaussianQuadrature(fPi, 0.0, 1.0, 3, truePi);
    assert(std::abs(resG3.value - truePi) < 1e-3);

    auto resG5 = IntegrationLab::gaussianQuadrature(fPi, 0.0, 1.0, 5, truePi);
    assert(std::abs(resG5.value - truePi) < 1e-5);

    // Known polynomial integral: int_0^2 x^3 dx = [x^4 / 4]_0^2 = 16 / 4 = 4.0
    // 2-point Gauss quadrature is exact for polynomials up to degree 2*2 - 1 = 3!
    auto fCubic = [](double x){ return x * x * x; };
    auto resG2Exact = IntegrationLab::gaussianQuadrature(fCubic, 0.0, 2.0, 2, 4.0);
    assert(std::abs(resG2Exact.value - 4.0) < 1e-12);

    std::cout << "  -> Unit II Numerical Integration tests PASSED." << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "    EES UNIT II AUTOMATED TEST SUITE    " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testInterpolation();
        testIntegration();

        std::cout << "\n========================================" << std::endl;
        std::cout << "    ALL UNIT II TESTS PASSED (2/2)      " << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "UNIT II TEST FAILURE: " << ex.what() << std::endl;
        return 1;
    }
}
