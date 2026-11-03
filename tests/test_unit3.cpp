#include "ees/unit3/LinearSystems.hpp"
#include "ees/unit3/Splines.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace ees::unit3;

void testLinearSystems() {
    std::cout << "[TEST] Running Unit III Linear Systems & Pivoting tests..." << std::endl;

    // System:
    //  2x +  y -  z =  8
    // -3x -  y + 2z = -11
    // -2x +  y + 2z = -3
    // Exact solution: x = 2, y = 3, z = -1
    ees::core::Matrix A{
        { 2.0,  1.0, -1.0},
        {-3.0, -1.0,  2.0},
        {-2.0,  1.0,  2.0}
    };
    ees::core::Vector b{8.0, -11.0, -3.0};

    // 1. Gauss Elimination with Partial Pivoting
    auto resGauss = LinearAlgebraLab::gaussElimination(A, b, true);
    assert(resGauss.success);
    assert(std::abs(resGauss.solution[0] - 2.0) < 1e-10);
    assert(std::abs(resGauss.solution[1] - 3.0) < 1e-10);
    assert(std::abs(resGauss.solution[2] - (-1.0)) < 1e-10);

    // 2. Gauss-Jordan Elimination
    auto resGJ = LinearAlgebraLab::gaussJordan(A, b);
    assert(resGJ.success);
    assert(std::abs(resGJ.solution[0] - 2.0) < 1e-10);
    assert(std::abs(resGJ.solution[1] - 3.0) < 1e-10);
    assert(std::abs(resGJ.solution[2] - (-1.0)) < 1e-10);

    // 3. Matrix Inverse: A * A^-1 == I
    auto invA = LinearAlgebraLab::inverse(A);
    auto I = A * invA;
    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 3; ++c) {
            double expected = (r == c) ? 1.0 : 0.0;
            assert(std::abs(I(r, c) - expected) < 1e-10);
        }
    }

    std::cout << "  -> Unit III Linear Systems tests PASSED." << std::endl;
}

void testMatrixFactorizations() {
    std::cout << "[TEST] Running Unit III Matrix Factorization tests..." << std::endl;

    ees::core::Matrix A{
        { 4.0,  3.0},
        { 6.0,  3.0}
    };

    // 1. Doolittle
    auto doo = LinearAlgebraLab::doolittleLU(A);
    assert(doo.success);
    assert(doo.reconstructionError < 1e-12);
    // Doolittle: L has 1s on diagonal
    assert(std::abs(doo.L(0, 0) - 1.0) < 1e-12);
    assert(std::abs(doo.L(1, 1) - 1.0) < 1e-12);

    // 2. Crout
    auto crout = LinearAlgebraLab::croutLU(A);
    assert(crout.success);
    assert(crout.reconstructionError < 1e-12);
    // Crout: U has 1s on diagonal
    assert(std::abs(crout.U(0, 0) - 1.0) < 1e-12);
    assert(std::abs(crout.U(1, 1) - 1.0) < 1e-12);

    // 3. Cholesky on Symmetric Positive Definite matrix
    // A = [ 4  12 ]
    //     [ 12 45 ] => L = [ 2 0 ], L * L^T = [ 4  12 ]
    //                      [ 6 3 ]            [ 12 45 ]
    ees::core::Matrix Aspd{
        { 4.0, 12.0},
        {12.0, 45.0}
    };
    auto chol = LinearAlgebraLab::cholesky(Aspd);
    assert(chol.success);
    assert(chol.reconstructionError < 1e-12);
    assert(std::abs(chol.L(0, 0) - 2.0) < 1e-12);
    assert(std::abs(chol.L(1, 0) - 6.0) < 1e-12);
    assert(std::abs(chol.L(1, 1) - 3.0) < 1e-12);

    std::cout << "  -> Unit III Factorization tests PASSED." << std::endl;
}

void testPowerMethod() {
    std::cout << "[TEST] Running Unit III Power Method tests..." << std::endl;

    // Matrix with eigenvalues 4 and 1
    // A = [ 2  2 ]
    //     [ 1  3 ]
    // Characteristic eq: (2-L)(3-L) - 2 = L^2 - 5L + 4 = (L - 4)(L - 1) = 0 => lambda_max = 4.0
    ees::core::Matrix A{
        {2.0, 2.0},
        {1.0, 3.0}
    };
    auto eig = LinearAlgebraLab::powerMethod(A, ees::core::Vector{1.0, 0.0}, 1e-7);
    assert(eig.success);
    assert(std::abs(eig.eigenvalue - 4.0) < 1e-5);

    std::cout << "  -> Unit III Power Method tests PASSED." << std::endl;
}

void testSplines() {
    std::cout << "[TEST] Running Unit III Spline tests..." << std::endl;

    ees::core::Vector x{0.0, 1.0, 2.0, 3.0};
    ees::core::Vector y{0.0, 1.0, 8.0, 27.0}; // y = x^3

    // Natural Cubic Spline must interpolate knots exactly: S(x_i) == y_i
    for (size_t i = 0; i < x.size(); ++i) {
        auto res = SplineLab::naturalCubicSpline(x, y, x[i]);
        assert(std::abs(res.interpolatedY - y[i]) < 1e-10);
    }

    // Evaluate between knots: x = 1.5
    auto resMid = SplineLab::naturalCubicSpline(x, y, 1.5);
    assert(resMid.interpolatedY > 1.0 && resMid.interpolatedY < 8.0);

    // Linear spline at x = 1.5 is midpoint of 1 and 8: 4.5
    auto linMid = SplineLab::linearSpline(x, y, 1.5);
    assert(std::abs(linMid.interpolatedY - 4.5) < 1e-10);

    // B-Spline evaluation check
    auto bsplineCurve = SplineLab::evaluateBSpline(x, y, 20);
    assert(bsplineCurve.size() == 20);

    std::cout << "  -> Unit III Spline tests PASSED." << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "   EES UNIT III AUTOMATED TEST SUITE    " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testLinearSystems();
        testMatrixFactorizations();
        testPowerMethod();
        testSplines();

        std::cout << "\n========================================" << std::endl;
        std::cout << "   ALL UNIT III TESTS PASSED (4/4)      " << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "UNIT III TEST FAILURE: " << ex.what() << std::endl;
        return 1;
    }
}

