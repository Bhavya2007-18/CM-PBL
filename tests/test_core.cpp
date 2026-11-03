#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include "ees/core/NumericalResult.hpp"
#include "ees/core/ErrorAnalysis.hpp"
#include "ees/core/ExpressionParser.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace ees::core;

void testVector() {
    std::cout << "[TEST] Running Vector tests..." << std::endl;

    Vector v1{1.0, 2.0, 3.0};
    Vector v2{4.0, 5.0, 6.0};

    // Addition & Subtraction
    Vector vAdd = v1 + v2;
    assert(std::abs(vAdd[0] - 5.0) < 1e-12);
    assert(std::abs(vAdd[1] - 7.0) < 1e-12);
    assert(std::abs(vAdd[2] - 9.0) < 1e-12);

    Vector vSub = v2 - v1;
    assert(std::abs(vSub[0] - 3.0) < 1e-12);
    assert(std::abs(vSub[1] - 3.0) < 1e-12);
    assert(std::abs(vSub[2] - 3.0) < 1e-12);

    // Dot product: 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    assert(std::abs(v1.dot(v2) - 32.0) < 1e-12);

    // Cross product: (1,2,3) x (4,5,6) = (-3, 6, -3)
    Vector vCross = v1.cross(v2);
    assert(std::abs(vCross[0] - (-3.0)) < 1e-12);
    assert(std::abs(vCross[1] - 6.0) < 1e-12);
    assert(std::abs(vCross[2] - (-3.0)) < 1e-12);

    // Norms
    Vector v3{3.0, -4.0};
    assert(std::abs(v3.norm1() - 7.0) < 1e-12);
    assert(std::abs(v3.norm2() - 5.0) < 1e-12);
    assert(std::abs(v3.normInf() - 4.0) < 1e-12);

    // Linspace
    Vector lin = Vector::linspace(0.0, 10.0, 5);
    assert(lin.size() == 5);
    assert(std::abs(lin[0] - 0.0) < 1e-12);
    assert(std::abs(lin[2] - 5.0) < 1e-12);
    assert(std::abs(lin[4] - 10.0) < 1e-12);

    std::cout << "  -> Vector tests PASSED." << std::endl;
}

void testMatrix() {
    std::cout << "[TEST] Running Matrix tests..." << std::endl;

    Matrix A{
        {1.0, 2.0},
        {3.0, 4.0}
    };

    Matrix B{
        {2.0, 0.0},
        {1.0, 2.0}
    };

    // Matrix multiplication:
    // [1 2] * [2 0] = [1*2+2*1  1*0+2*2] = [4 4]
    // [3 4]   [1 2]   [3*2+4*1  3*0+4*2]   [10 8]
    Matrix C = A * B;
    assert(std::abs(C(0, 0) - 4.0) < 1e-12);
    assert(std::abs(C(0, 1) - 4.0) < 1e-12);
    assert(std::abs(C(1, 0) - 10.0) < 1e-12);
    assert(std::abs(C(1, 1) - 8.0) < 1e-12);

    // Matrix * Vector
    Vector x{1.0, 2.0};
    Vector Ax = A * x; // [1*1+2*2, 3*1+4*2] = [5, 11]
    assert(std::abs(Ax[0] - 5.0) < 1e-12);
    assert(std::abs(Ax[1] - 11.0) < 1e-12);

    // Transpose
    Matrix At = A.transpose();
    assert(std::abs(At(0, 1) - 3.0) < 1e-12);
    assert(std::abs(At(1, 0) - 2.0) < 1e-12);

    // Determinant: 1*4 - 2*3 = -2
    assert(std::abs(A.determinant() - (-2.0)) < 1e-12);

    // Row operations
    Matrix M{
        {2.0, 4.0},
        {1.0, 3.0}
    };
    M.swapRows(0, 1);
    assert(std::abs(M(0, 0) - 1.0) < 1e-12);
    assert(std::abs(M(1, 0) - 2.0) < 1e-12);

    M.addRowMultiple(1, 0, -2.0); // R1 = R1 - 2*R0 => [2, 4] - 2*[1, 3] = [0, -2]
    assert(std::abs(M(1, 0) - 0.0) < 1e-12);
    assert(std::abs(M(1, 1) - (-2.0)) < 1e-12);

    std::cout << "  -> Matrix tests PASSED." << std::endl;
}

void testErrorAnalysis() {
    std::cout << "[TEST] Running Error Analysis tests..." << std::endl;

    double trueVal = 10.0 / 3.0; // 3.3333333333...
    double approxVal = 3.333;

    double absErr = ErrorAnalysis::absoluteError(trueVal, approxVal);
    double relErr = ErrorAnalysis::relativeError(trueVal, approxVal);
    double pctErr = ErrorAnalysis::percentageError(trueVal, approxVal);

    assert(std::abs(absErr - 0.0003333333333333) < 1e-10);
    assert(std::abs(relErr - 0.0001) < 1e-6);
    assert(std::abs(pctErr - 0.01) < 1e-4);

    int sigDigits = ErrorAnalysis::significantDigits(relErr);
    assert(sigDigits >= 3 && sigDigits <= 4);

    // IEEE 754 inspection
    auto info = ErrorAnalysis::inspectDouble(1.0);
    assert(info.sign == 0);
    assert(info.exponent == 1023); // Unbiased 0
    assert(info.mantissa == 0);

    // Catastrophic cancellation demo
    double largeX = 1e12;
    double direct = ErrorAnalysis::cancellationDirect(largeX);
    double stable = ErrorAnalysis::cancellationStable(largeX);
    assert(std::abs(stable - 5e-7) < 1e-12);
    assert(std::abs(direct - stable) < 1e-10);

    std::cout << "  -> Error Analysis tests PASSED." << std::endl;
}

void testExpressionParser() {
    std::cout << "[TEST] Running Expression Parser tests..." << std::endl;

    ExpressionParser p1("x^3 - 4*x - 9");
    assert(p1.isValid());
    // For x = 3: 27 - 12 - 9 = 6
    assert(std::abs(p1.evaluate(3.0) - 6.0) < 1e-12);

    // Derivative: d/dx (x^3 - 4x - 9) = 3x^2 - 4. At x = 3: 3*(9) - 4 = 23
    double dCentral = p1.derivativeCentral(3.0);
    assert(std::abs(dCentral - 23.0) < 1e-4);

    // 2nd derivative: 6x. At x = 3: 18
    double d2 = p1.secondDerivative(3.0);
    assert(std::abs(d2 - 18.0) < 1e-3);

    // Trig and transcendental expressions
    ExpressionParser p2("sin(x) + cos(x)");
    assert(std::abs(p2.evaluate(0.0) - 1.0) < 1e-12);

    ExpressionParser p3("exp(x) - 1");
    assert(std::abs(p3.evaluate(0.0) - 0.0) < 1e-12);

    // Multivariate expression
    ExpressionParser pMulti("x^2 + y^2");
    assert(std::abs(pMulti.evaluate(3.0, 4.0) - 25.0) < 1e-12);

    // Multivariate gradient at (3, 4) of x^2 + y^2 is [2x, 2y] = [6, 8]
    Vector pt{3.0, 4.0};
    Vector grad = pMulti.gradient(pt, {"x", "y"});
    assert(std::abs(grad[0] - 6.0) < 1e-4);
    assert(std::abs(grad[1] - 8.0) < 1e-4);

    std::cout << "  -> Expression Parser tests PASSED." << std::endl;
}

void testNumericalResult() {
    std::cout << "[TEST] Running Numerical Result tests..." << std::endl;

    NumericalResult res;
    res.success = true;
    res.solution = 2.70653;
    res.status = ConvergenceStatus::Converged;
    res.convergenceRate = ConvergenceRate::Quadratic;
    res.convergenceRateExplanation = "Newton-Raphson exhibits quadratic asymptotic error reduction.";
    res.executionTimeMs = 0.42;

    IterationRecord r1{1, 2.0, Vector(), -5.0, 0.70653, 0.261, {{"f'(x)", 8.0}}, "Initial step"};
    IterationRecord r2{2, 2.70653, Vector(), 0.001, 0.0001, 0.00003, {{"f'(x)", 18.2}}, "Converged"};
    res.addRecord(r1);
    res.addRecord(r2);

    std::string summary = res.summary();
    assert(summary.find("Converged") != std::string::npos);
    assert(summary.find("Quadratic") != std::string::npos);

    std::string table = res.tableString({"f'(x)"});
    assert(table.find("Iter") != std::string::npos);
    assert(table.find("Initial step") != std::string::npos);

    std::cout << "  -> Numerical Result tests PASSED." << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "   EES CORE MATHEMATICAL TESTS SUITE    " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testVector();
        testMatrix();
        testErrorAnalysis();
        testExpressionParser();
        testNumericalResult();

        std::cout << "\n========================================" << std::endl;
        std::cout << "   ALL CORE UNIT TESTS PASSED (5/5)     " << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "TEST FAILURE: " << ex.what() << std::endl;
        return 1;
    }
}

