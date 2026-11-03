#include "ees/ui/Workbench.hpp"
#include "ees/ui/Visualizer.hpp"
#include "ees/unit1/Foundations.hpp"
#include "ees/unit1/RootFinding.hpp"
#include "ees/unit1/Optimization.hpp"
#include "ees/unit2/Interpolation.hpp"
#include "ees/unit2/Integration.hpp"
#include "ees/unit3/LinearSystems.hpp"
#include "ees/unit3/Splines.hpp"
#include "ees/unit4/ODESolvers.hpp"
#include "ees/unit4/PDESolvers.hpp"
#include "ees/comparison/MethodComparator.hpp"
#include "ees/problems/ProblemLibrary.hpp"
#include "ees/advisor/MethodAdvisor.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

namespace ees::ui {

Workbench::Workbench() = default;

void Workbench::printHeader(const std::string& title) {
    std::cout << "\n\033[1;36m========================================================================================\033[0m\n";
    std::cout << "\033[1;37m   EES — ENGINEERING EQUATION SOLVER\033[0m\n";
    std::cout << "\033[1;33m   " << title << "\033[0m\n";
    std::cout << "\033[1;36m========================================================================================\033[0m\n";
}

void Workbench::pause() {
    std::cout << "\n\033[90mPress Enter to return to menu...\033[0m";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int Workbench::readInt(int minVal, int maxVal) {
    int choice = -1;
    while (true) {
        std::cout << "\033[1;32mSelect option [" << minVal << "-" << maxVal << "]: \033[0m";
        if (std::cin >> choice && choice >= minVal && choice <= maxVal) {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return choice;
        }
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "\033[1;31mInvalid input. Please enter a valid number.\033[0m\n";
    }
}

double Workbench::readDouble(const std::string& prompt, double defaultVal) {
    std::cout << prompt << " [" << defaultVal << "]: ";
    std::string input;
    std::getline(std::cin, input);
    if (input.empty()) return defaultVal;
    try {
        return std::stod(input);
    } catch (...) {
        return defaultVal;
    }
}

std::string Workbench::readString(const std::string& prompt, const std::string& defaultVal) {
    std::cout << prompt << " [" << defaultVal << "]: ";
    std::string input;
    std::getline(std::cin, input);
    if (input.empty()) return defaultVal;
    return input;
}

void Workbench::run() {
    while (true) {
        printHeader("MAIN WORKBENCH DASHBOARD");
        std::cout << "\033[1m  [1]\033[0m UNIT I:  Foundations, Error Analysis & Computer Arithmetic Lab\n";
        std::cout << "\033[1m  [2]\033[0m UNIT I:  Root-Finding Laboratory (Bisection, Newton, Secant)\n";
        std::cout << "\033[1m  [3]\033[0m UNIT I:  Optimization (Fibonacci, Golden Section, Nelder-Mead)\n";
        std::cout << "\033[1m  [4]\033[0m UNIT II: Interpolation & Finite Differences Laboratory\n";
        std::cout << "\033[1m  [5]\033[0m UNIT II: Numerical Integration & Quadrature Laboratory\n";
        std::cout << "\033[1m  [6]\033[0m UNIT III: Linear Systems, Pivoting, LU & Cholesky Decomposition\n";
        std::cout << "\033[1m  [7]\033[0m UNIT III: Dominant Eigenvalues (Power Method)\n";
        std::cout << "\033[1m  [8]\033[0m UNIT III: Spline Approximation (Linear, Quadratic, Natural Cubic, B-Spline)\n";
        std::cout << "\033[1m  [9]\033[0m UNIT IV: Numerical Solution of ODEs (Euler, RK4, Predictor-Corrector)\n";
        std::cout << "\033[1m [10]\033[0m UNIT IV: Numerical Solution of PDEs (Heat, Wave, Laplace/Poisson)\n";
        std::cout << "\033[1m [11]\033[0m METHOD COMPARISON BENCHMARK ENGINE\n";
        std::cout << "\033[1m [12]\033[0m ENGINEERING PROBLEM CATALOG (Mechanical, Civil, Electrical, Thermal)\n";
        std::cout << "\033[1m [13]\033[0m NUMERICAL METHOD ADVISOR (Intelligent Recommendation Engine)\n";
        std::cout << "\033[1;32m [14]\033[0m ALL-IN-ONE PBL VIVA DEMONSTRATION MODE (One-Click Showcase)\033[0m\n";
        std::cout << "\033[1m  [0]\033[0m Exit EES\n\n";

        int opt = readInt(0, 14);
        if (opt == 0) {
            std::cout << "\n\033[1;33mExiting EES Engineering Equation Solver. Good luck with your PBL presentation!\033[0m\n";
            break;
        }

        switch (opt) {
            case 1:  runUnit1Foundations(); break;
            case 2:  runUnit1Roots(); break;
            case 3:  runUnit1Optimization(); break;
            case 4:  runUnit2Interpolation(); break;
            case 5:  runUnit2Integration(); break;
            case 6:  runUnit3LinearSystems(); break;
            case 7:  runUnit3Eigenvalues(); break;
            case 8:  runUnit3Splines(); break;
            case 9:  runUnit4ODE(); break;
            case 10: runUnit4PDE(); break;
            case 11: runMethodComparison(); break;
            case 12: runProblemLibrary(); break;
            case 13: runMethodAdvisor(); break;
            case 14: runDemoMode(); break;
            default: break;
        }
    }
}

void Workbench::runUnit1Foundations() {
    printHeader("UNIT I — NUMERICAL FOUNDATIONS & COMPUTER ARITHMETIC");
    std::cout << "1. Computer Arithmetic & Catastrophic Cancellation Lab\n";
    std::cout << "2. Taylor Series Polynomial Expansion\n";
    std::cout << "3. Rolle's Theorem & Mean Value Theorem Verification\n";
    std::cout << "0. Back to Main Menu\n\n";

    int choice = readInt(0, 3);
    if (choice == 1) {
        std::cout << unit1::FoundationsLab::runArithmeticStabilityLab() << "\n";
    } else if (choice == 2) {
        std::string expr = readString("Enter function f(x)", "exp(x)");
        double x0 = readDouble("Expansion point x0", 0.0);
        double x = readDouble("Evaluation point x", 1.0);
        int maxOrder = static_cast<int>(readDouble("Order / Number of terms", 5));
        double trueVal = readDouble("Analytical True Value (for error tracking)", std::exp(1.0));

        auto res = unit1::FoundationsLab::computeTaylorSeriesExpr(expr, x0, x, maxOrder, trueVal);
        std::cout << "\n" << res.summary() << "\n";
    } else if (choice == 3) {
        auto f = [](double x){ return x * (x - 2.0); };
        std::cout << "\nRolle's Theorem Demo on f(x) = x*(x - 2) on [0, 2]:\n";
        auto rRolle = unit1::FoundationsLab::verifyRollesTheorem(f, 0.0, 2.0);
        std::cout << rRolle.explanation << "\n";

        auto fCubic = [](double x){ return x * x * x; };
        std::cout << "\nMean Value Theorem Demo on f(x) = x^3 on [0, 3]:\n";
        auto rMVT = unit1::FoundationsLab::verifyMeanValueTheorem(fCubic, 0.0, 3.0);
        std::cout << rMVT.explanation << "\n";
    }
    pause();
}

void Workbench::runUnit1Roots() {
    printHeader("UNIT I — ROOT-FINDING LABORATORY");
    std::string expr = readString("Enter nonlinear function f(x)", "x^3 - 4*x - 9");
    double a = readDouble("Lower bracket / guess 1 (a)", 2.0);
    double b = readDouble("Upper bracket / guess 2 (b)", 3.0);
    double tol = readDouble("Convergence Tolerance", 1e-7);

    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };

    std::cout << "\nSolving via Bisection, Newton-Raphson, and Secant...\n";
    auto rBis = unit1::RootFindingLab::bisection(f, a, b, tol);
    auto rNew = unit1::RootFindingLab::newtonRaphson(f, 0.5 * (a + b), tol);
    auto rSec = unit1::RootFindingLab::secant(f, a, b, tol);

    std::cout << "\n" << rBis.summary() << "\n";
    std::cout << rBis.tableString({"a", "b", "interval_width"}) << "\n";

    std::cout << "\n" << rNew.summary() << "\n";
    std::cout << rNew.tableString({"f'(x)"}) << "\n";

    // Visualizations
    std::cout << Visualizer::plotFunction(f, a - 0.5, b + 0.5, 65, 15, {rNew.solution}, "Root Location on f(x)") << "\n";
    std::cout << Visualizer::plotConvergence(rNew.history, 65, 14, "Newton-Raphson Quadratic Convergence") << "\n";
    std::cout << Visualizer::plotIntervalShrinkage(rBis.history, 50) << "\n";

    std::cout << "\n" << unit1::RootFindingLab::compareRootMethods(f, a, b, 0.5 * (a + b), tol) << "\n";
    pause();
}

void Workbench::runUnit1Optimization() {
    printHeader("UNIT I — UNCONSTRAINED OPTIMIZATION LABORATORY");
    std::cout << "1. 1D Optimization (Fibonacci vs Golden Section vs Newton 1D)\n";
    std::cout << "2. Multivariate Optimization (Steepest Descent vs Nelder-Mead Simplex)\n";
    int opt = readInt(1, 2);

    if (opt == 1) {
        auto f1 = [](double x){ return (x - 3.0) * (x - 3.0) + 4.0; };
        std::cout << "Target Function: f(x) = (x - 3)^2 + 4, Minimum at x = 3.0, f(3) = 4.0\n";
        double a = readDouble("Interval lower bound a", 0.0);
        double b = readDouble("Interval upper bound b", 5.0);

        auto rFib = unit1::OptimizationLab::fibonacciSearch(f1, a, b, 1e-4);
        auto rGold = unit1::OptimizationLab::goldenSectionSearch(f1, a, b, 1e-4);
        auto rNewt = unit1::OptimizationLab::newtonsMethod1D(f1, a, 1e-6);

        std::cout << "\n" << rGold.summary() << "\n";
        std::cout << rGold.tableString({"a", "b", "width"}) << "\n";
        std::cout << Visualizer::plotIntervalShrinkage(rGold.history, 50) << "\n";
    } else {
        auto fMulti = [](const core::Vector& v){
            return (v[0] - 2.0) * (v[0] - 2.0) + (v[1] - 3.0) * (v[1] - 3.0);
        };
        std::cout << "Target Function: f(x, y) = (x - 2)^2 + (y - 3)^2, Global Minimum at (2, 3)\n";
        core::Vector start{0.0, 0.0};

        auto rSD = unit1::OptimizationLab::steepestDescent(fMulti, start, 0.2, 1e-4);
        auto rNM = unit1::OptimizationLab::nelderMead(fMulti, start, 0.5, 1e-4);

        std::cout << "\nSteepest Descent Result:\n" << rSD.summary() << "\n";
        std::cout << "Nelder-Mead Simplex Result:\n" << rNM.summary() << "\n";
    }
    pause();
}

void Workbench::runUnit2Interpolation() {
    printHeader("UNIT II — INTERPOLATION & FINITE DIFFERENCES LABORATORY");
    core::Vector x{1.0, 2.0, 3.0, 4.0, 5.0};
    core::Vector y{1.0, 8.0, 27.0, 64.0, 125.0}; // y = x^3

    std::cout << "Dataset: x = [1, 2, 3, 4, 5], y = [1, 8, 27, 64, 125] (Cubic Data)\n";
    auto diffTbl = unit2::InterpolationLab::computeForwardDifferenceTable(x, y);
    std::cout << "\n" << diffTbl.toString() << "\n";

    double evalX = readDouble("Enter evaluation point x", 2.5);
    std::cout << "\n" << unit2::InterpolationLab::compareInterpolation(x, y, evalX) << "\n";

    auto lag = unit2::InterpolationLab::lagrange(x, y, evalX);
    std::cout << lag.summary() << "\n";
    pause();
}

void Workbench::runUnit2Integration() {
    printHeader("UNIT II — NUMERICAL INTEGRATION & QUADRATURE LABORATORY");
    auto fPi = [](double x){ return 4.0 / (1.0 + x * x); };
    const double truePi = 3.14159265358979323846;

    std::cout << "Standard Integral: int_0^1 (4 / (1 + x^2)) dx = pi ≈ 3.1415926535...\n";
    std::cout << "\n" << unit2::IntegrationLab::compareIntegrationMethods(fPi, 0.0, 1.0, truePi) << "\n";

    auto rRomb = unit2::IntegrationLab::romberg(fPi, 0.0, 1.0, 5, 1e-9, truePi);
    std::cout << rRomb.summary() << "\n";
    pause();
}

void Workbench::runUnit3LinearSystems() {
    printHeader("UNIT III — LINEAR SYSTEMS & MATRIX FACTORIZATION");
    std::cout << "1. Solve Linear System via Gauss Elimination with Partial Pivoting\n";
    std::cout << "2. Solve via Gauss-Jordan (RREF & Matrix Inverse)\n";
    std::cout << "3. Matrix Factorization (Doolittle, Crout, Cholesky)\n";
    std::cout << "4. Numerical Pivoting & Stability Demonstration\n";
    int opt = readInt(1, 4);

    if (opt == 1 || opt == 2) {
        core::Matrix A{
            { 2.0,  1.0, -1.0},
            {-3.0, -1.0,  2.0},
            {-2.0,  1.0,  2.0}
        };
        core::Vector b{8.0, -11.0, -3.0};
        std::cout << "Matrix A:\n" << A.toString() << "\nRHS b:\n" << b.toString() << "\n\n";

        if (opt == 1) {
            auto sol = unit3::LinearAlgebraLab::gaussElimination(A, b, true);
            std::cout << sol.summary() << "\n";
            std::cout << "Row Operation Step History:\n";
            for (const auto& step : sol.stepHistory) {
                std::cout << "  Step " << step.stepNumber << ": " << step.operation << "\n";
            }
        } else {
            auto sol = unit3::LinearAlgebraLab::gaussJordan(A, b);
            std::cout << sol.summary() << "\n";
            auto inv = unit3::LinearAlgebraLab::inverse(A);
            std::cout << "Computed Matrix Inverse A^-1:\n" << inv.toString() << "\n";
        }
    } else if (opt == 3) {
        core::Matrix Aspd{
            { 4.0, 12.0},
            {12.0, 45.0}
        };
        std::cout << "Symmetric Positive-Definite Matrix A:\n" << Aspd.toString() << "\n\n";
        auto doo = unit3::LinearAlgebraLab::doolittleLU(Aspd);
        auto chol = unit3::LinearAlgebraLab::cholesky(Aspd);
        std::cout << doo.summary() << "\n\n" << chol.summary() << "\n";
    } else {
        std::cout << unit3::LinearAlgebraLab::runPivotingDemonstration() << "\n";
    }
    pause();
}

void Workbench::runUnit3Eigenvalues() {
    printHeader("UNIT III — POWER METHOD DOMINANT EIGENVALUE SOLVER");
    core::Matrix A{
        {2.0, 2.0},
        {1.0, 3.0}
    };
    std::cout << "Input Matrix A:\n" << A.toString() << "\n";
    auto eig = unit3::LinearAlgebraLab::powerMethod(A, core::Vector{1.0, 0.0}, 1e-7);
    std::cout << "\n" << eig.summary() << "\n";
    std::cout << Visualizer::plotConvergence(eig.history, 65, 14, "Power Method Rayleigh Quotient Convergence") << "\n";
    pause();
}

void Workbench::runUnit3Splines() {
    printHeader("UNIT III — SPLINE APPROXIMATION LABORATORY");
    core::Vector x{0.0, 1.0, 2.0, 3.0, 4.0};
    core::Vector y{0.0, 1.0, 8.0, 27.0, 64.0}; // y = x^3

    std::cout << "Knot Points: x = [0, 1, 2, 3, 4], y = [0, 1, 8, 27, 64]\n\n";
    double evalX = readDouble("Evaluation point x", 2.5);

    auto cubic = unit3::SplineLab::naturalCubicSpline(x, y, evalX);
    std::cout << "\n" << cubic.summary() << "\n";
    std::cout << unit3::SplineLab::compareSplines(x, y, evalX) << "\n";
    pause();
}

void Workbench::runUnit4ODE() {
    printHeader("UNIT IV — NUMERICAL SOLUTION OF ODES LABORATORY");
    auto f = [](double x, double y){ (void)x; return y; };
    auto exact = [](double x){ return std::exp(x); };

    std::cout << "IVP: dy/dx = y, y(0) = 1.0 on interval [0, 1]. Exact Solution: y(x) = exp(x)\n\n";
    double h = readDouble("Step size h", 0.1);

    std::cout << "\n" << unit4::ODESolverLab::compareODESolvers(f, 0.0, 1.0, 1.0, h, exact) << "\n";
    auto rRK4 = unit4::ODESolverLab::rungeKutta4(f, 0.0, 1.0, 1.0, h, exact);
    std::cout << rRK4.tableString(10) << "\n";
    pause();
}

void Workbench::runUnit4PDE() {
    printHeader("UNIT IV — NUMERICAL SOLUTION OF PDES LABORATORY");
    std::cout << "1. Parabolic: 1D Heat Equation (FTCS Scheme with CFL Stability)\n";
    std::cout << "2. Elliptic: 2D Laplace Equation (Gauss-Seidel with SOR)\n";
    std::cout << "3. Hyperbolic: 1D Wave Equation (Explicit Central Difference Scheme)\n";
    int opt = readInt(1, 3);

    if (opt == 1) {
        auto u0 = [](double x){ return 100.0 * std::sin(3.141592653589793 * x); };
        auto heatSol = unit4::PDESolverLab::solveHeatEquation1D(0.05, 1.0, 0.5, 21, 51, u0, 0.0, 0.0);
        std::cout << "\n" << heatSol.summary() << "\n";
        std::cout << heatSol.asciiHeatmap(50, 15) << "\n";
    } else if (opt == 2) {
        auto laplaceSol = unit4::PDESolverLab::solveLaplaceEquation2D(1.0, 1.0, 15, 15, 100.0, 0.0, 0.0, 0.0, 1.3, 1e-4, 500);
        std::cout << "\n" << laplaceSol.summary() << "\n";
        std::cout << laplaceSol.asciiHeatmap(45, 15) << "\n";
    } else {
        auto u0 = [](double x){ return std::sin(3.141592653589793 * x); };
        auto waveSol = unit4::PDESolverLab::solveWaveEquation1D(1.0, 1.0, 1.0, 21, 51, u0, nullptr, 0.0, 0.0);
        std::cout << "\n" << waveSol.summary() << "\n";
        std::cout << waveSol.asciiHeatmap(50, 15) << "\n";
    }
    pause();
}

void Workbench::runProblemLibrary() {
    printHeader("ENGINEERING PROBLEM LIBRARY");
    auto problems = problems::ProblemLibrary::getAllProblems();

    for (size_t i = 0; i < problems.size(); ++i) {
        std::cout << "[" << i + 1 << "] " << problems[i].id << ": " << problems[i].title << "\n";
    }
    std::cout << "[0] Return\n\n";

    int choice = readInt(0, static_cast<int>(problems.size()));
    if (choice == 0) return;

    const auto& p = problems[choice - 1];
    std::cout << "\n=== " << p.title << " (" << p.id << ") ===\n";
    std::cout << "Description:      " << p.description << "\n";
    std::cout << "Math Model:       " << p.mathematicalModel << "\n";
    std::cout << "Known Parameters: " << p.parametersExplanation << "\n";
    std::cout << "Recommended:      " << p.recommendedMethods << "\n";
    std::cout << "Benchmark Result: " << p.expectedResult << "\n\n";

    std::cout << "Executing C++ Numerical Methods Engine on this problem...\n";
    if (p.id == "MECH-01") {
        double v = problems::ProblemLibrary::solveVanDerWaals(5.066e6, 320.0, 8.314, 0.3643, 4.267e-5);
        std::cout << "\033[1;32mComputed Molar Volume: " << std::setprecision(8) << v << " m^3/mol\033[0m\n";
    } else if (p.id == "CIVIL-01") {
        auto F = problems::ProblemLibrary::solveTrussNodalForces();
        std::cout << "\033[1;32mComputed Member Forces F:\n" << F.toString(4) << " kN\033[0m\n";
    } else if (p.id == "ELEC-01") {
        double vd = problems::ProblemLibrary::solveDiodeOperatingPoint(5.0, 1000.0, 1e-12, 0.026);
        std::cout << "\033[1;32mComputed Diode Forward Voltage: " << vd << " V\033[0m\n";
    } else if (p.id == "GEN-01") {
        double mu = problems::ProblemLibrary::interpolateOilViscosity(50.0);
        std::cout << "\033[1;32mComputed Oil Viscosity at 50 C: " << mu << " mPa.s\033[0m\n";
    }
    pause();
}

void Workbench::runMethodComparison() {
    printHeader("METHOD COMPARISON BENCHMARK ENGINE");
    std::cout << "1. Root Finding (Bisection vs Newton-Raphson vs Secant)\n";
    std::cout << "2. Optimization (Fibonacci vs Golden Section vs Newton 1D)\n";
    std::cout << "3. Numerical Integration (Trapezoidal vs Simpson vs Romberg vs Gauss-Legendre)\n";
    std::cout << "4. ODE Solvers (Euler vs Heun vs RK4 vs Adams-Bashforth vs Milne)\n";
    int opt = readInt(1, 4);

    if (opt == 1) {
        auto f = [](double x){ return x*x*x - 4.0*x - 9.0; };
        auto rep = comparison::MethodComparator::compareRoots(f, 2.0, 3.0, 2.5, 1e-7, "x^3 - 4x - 9 = 0");
        std::cout << "\n" << rep.summary() << "\n";
    } else if (opt == 2) {
        auto f = [](double x){ return (x - 3.0)*(x - 3.0) + 4.0; };
        auto rep = comparison::MethodComparator::compare1DOptimization(f, 0.0, 5.0, 1.0, 1e-5, "(x - 3)^2 + 4");
        std::cout << "\n" << rep.summary() << "\n";
    } else if (opt == 3) {
        auto f = [](double x){ return 4.0 / (1.0 + x*x); };
        auto rep = comparison::MethodComparator::compareIntegration(f, 0.0, 1.0, 3.141592653589793, "4 / (1 + x^2)");
        std::cout << "\n" << rep.summary() << "\n";
    } else {
        auto f = [](double x, double y){ (void)x; return y; };
        auto rep = comparison::MethodComparator::compareODESolvers(f, 0.0, 1.0, 1.0, 0.1, [](double x){ return std::exp(x); }, "dy/dx = y");
        std::cout << "\n" << rep.summary() << "\n";
    }
    pause();
}

void Workbench::runMethodAdvisor() {
    printHeader("NUMERICAL METHOD ADVISOR");
    std::cout << "Specify your problem characteristics:\n";
    std::cout << "Problem Category:\n";
    std::cout << "1. Root Finding\n2. Linear Systems\n3. Integration\n4. Interpolation\n5. ODE\n";
    int cat = readInt(1, 5);

    advisor::ProblemProfile prof;
    if (cat == 1) {
        prof.category = advisor::ProblemCategory::RootFinding;
        std::cout << "Is a bracketing interval [a, b] with sign change known? (1: Yes, 0: No): ";
        prof.hasBracketingInterval = (readInt(0, 1) == 1);
        std::cout << "Is the analytical derivative f'(x) computable? (1: Yes, 0: No): ";
        prof.hasAnalyticDerivative = (readInt(0, 1) == 1);
    } else if (cat == 2) {
        prof.category = advisor::ProblemCategory::LinearSystem;
        std::cout << "Is matrix A Symmetric and Positive Definite? (1: Yes, 0: No): ";
        prof.isMatrixSymmetricPositiveDefinite = (readInt(0, 1) == 1);
    } else if (cat == 3) {
        prof.category = advisor::ProblemCategory::NumericalIntegration;
        std::cout << "Is the integrand smooth and continuous? (1: Yes, 0: No): ";
        prof.isSmoothContinuous = (readInt(0, 1) == 1);
    } else if (cat == 4) {
        prof.category = advisor::ProblemCategory::Interpolation;
        std::cout << "Are the data points equispaced? (1: Yes, 0: No): ";
        prof.isEquispacedData = (readInt(0, 1) == 1);
    } else {
        prof.category = advisor::ProblemCategory::ODEInitialValue;
    }

    auto rec = advisor::MethodAdvisor::advise(prof);
    std::cout << "\n" << rec.toString() << "\n";
    pause();
}

void Workbench::runDemoMode(bool pauseAtEnd) {
    printHeader("ALL-IN-ONE PBL VIVA DEMONSTRATION SHOWCASE");
    std::cout << "\033[1;32mExecuting comprehensive demonstration across all 9 syllabus milestones...\033[0m\n\n";

    // 1. Root Finding
    std::cout << "\033[1;33m[1/9] ROOT FINDING (Newton-Raphson vs Bisection):\033[0m\n";
    auto fRoot = [](double x){ return x*x*x - 4.0*x - 9.0; };
    auto rNew = unit1::RootFindingLab::newtonRaphson(fRoot, 2.5, 1e-7);
    std::cout << "  Root of x^3 - 4x - 9 = 0: " << std::setprecision(8) << rNew.solution
              << " in " << rNew.iterations << " iterations (Order: Quadratic)\n";

    // 2. Optimization
    std::cout << "\033[1;33m[2/9] MULTIVARIATE OPTIMIZATION (Nelder-Mead Simplex):\033[0m\n";
    auto fOpt = [](const core::Vector& v){ return (v[0]-2.0)*(v[0]-2.0) + (v[1]-3.0)*(v[1]-3.0); };
    auto rNM = unit1::OptimizationLab::nelderMead(fOpt, core::Vector{0.0, 0.0}, 0.5, 1e-4);
    std::cout << "  Minimum found at: " << rNM.solutionVector.toString(4) << " in " << rNM.iterations << " cycles\n";

    // 3. Interpolation
    std::cout << "\033[1;33m[3/9] INTERPOLATION (Gregory-Newton Forward):\033[0m\n";
    core::Vector xi{1.0, 2.0, 3.0, 4.0};
    core::Vector yi{0.0, 5.0, 22.0, 57.0};
    auto rGN = unit2::InterpolationLab::gregoryNewtonForward(xi, yi, 2.5);
    std::cout << "  Interpolated y(2.5): " << rGN.interpolatedY << " (Exact: 11.625)\n";

    // 4. Numerical Integration
    std::cout << "\033[1;33m[4/9] NUMERICAL INTEGRATION (Gauss-Legendre 5-Point vs Romberg):\033[0m\n";
    auto fInt = [](double x){ return 4.0 / (1.0 + x*x); };
    auto rG5 = unit2::IntegrationLab::gaussianQuadrature(fInt, 0.0, 1.0, 5, 3.141592653589793);
    std::cout << "  Computed pi integral: " << std::setprecision(10) << rG5.value
              << " (Abs Error: " << std::scientific << rG5.absError << ")\n";

    // 5. Linear System
    std::cout << "\033[1;33m[5/9] LINEAR SYSTEM (Gauss Elimination with Partial Pivoting):\033[0m\n";
    core::Matrix A{{2.0, 1.0, -1.0}, {-3.0, -1.0, 2.0}, {-2.0, 1.0, 2.0}};
    core::Vector b{8.0, -11.0, -3.0};
    auto rGauss = unit3::LinearAlgebraLab::gaussElimination(A, b, true);
    std::cout << "  Solution x: " << rGauss.solution.toString(4) << " (Det: " << rGauss.determinant << ")\n";

    // 6. Dominant Eigenvalue
    std::cout << "\033[1;33m[6/9] EIGENVALUES (Power Method):\033[0m\n";
    core::Matrix Aeig{{2.0, 2.0}, {1.0, 3.0}};
    auto rEig = unit3::LinearAlgebraLab::powerMethod(Aeig, core::Vector{1.0, 0.0});
    std::cout << "  Dominant Eigenvalue: " << std::fixed << std::setprecision(6) << rEig.eigenvalue
              << " (Vector: " << rEig.eigenvector.toString(4) << ")\n";

    // 7. Natural Cubic Spline
    std::cout << "\033[1;33m[7/9] NATURAL CUBIC SPLINES (C2 Smoothness):\033[0m\n";
    core::Vector xs{0.0, 1.0, 2.0, 3.0};
    core::Vector ys{0.0, 1.0, 8.0, 27.0};
    auto rSpline = unit3::SplineLab::naturalCubicSpline(xs, ys, 1.5);
    std::cout << "  Spline interpolation at knot midpoint 1.5: " << rSpline.interpolatedY << "\n";

    // 8. ODE Solver
    std::cout << "\033[1;33m[8/9] ODE SOLUTION (Classic RK4):\033[0m\n";
    auto rRK4 = unit4::ODESolverLab::rungeKutta4([](double x, double y){ (void)x; return y; }, 0.0, 1.0, 1.0, 0.1, [](double x){ return std::exp(x); });
    std::cout << "  RK4 Final Value: " << rRK4.y.toStdVector().back() << " (Global Error: " << std::scientific << rRK4.maxError << ")\n";

    // 9. PDE Solver
    std::cout << "\033[1;33m[9/9] PDE SOLUTION (1D Heat Equation Diffusion):\033[0m\n";
    auto rHeat = unit4::PDESolverLab::solveHeatEquation1D(0.1, 1.0, 0.5, 11, 21, [](double x){ return std::sin(3.141592653589793*x); });
    std::cout << "  Grid Generated: " << rHeat.grid.rows() << " x " << rHeat.grid.cols() << " [CFL r = " << rHeat.parameterR << "]\n";
    std::cout << rHeat.asciiHeatmap(40, 10) << "\n";

    std::cout << "\033[1;32mAll 9 syllabus milestones verified and demonstrated successfully.\033[0m\n";
    if (pauseAtEnd) {
        pause();
    }
}

} // namespace ees::ui
