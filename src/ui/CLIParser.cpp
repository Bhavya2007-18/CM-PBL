#include "ees/ui/CLIParser.hpp"
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
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>

namespace ees::ui {

void CLIApp::printVersion() {
    std::cout << "EES (Engineering Equation Solver) v1.0.0\n";
    std::cout << "Numerical Methods & Engineering Analysis Workbench\n";
    std::cout << "Pure C++17 Core Engine — Zero External Dependencies\n";
}

void CLIApp::printUsage() {
    printVersion();
    std::cout << "\n\033[1mUSAGE:\033[0m\n";
    std::cout << "  ees [command] [options]\n";
    std::cout << "  ees interactive          Launch the interactive scientific workbench\n";
    std::cout << "  ees demo                 Run the complete 9-milestone PBL Viva Showcase\n\n";

    std::cout << "\033[1mAVAILABLE COMMANDS:\033[0m\n";
    std::cout << "  \033[1;36mroot\033[0m          Solve nonlinear equation f(x) = 0\n";
    std::cout << "  \033[1;36moptimize\033[0m      Find minimum/maximum of function f(x)\n";
    std::cout << "  \033[1;36mintegrate\033[0m     Compute definite integral of f(x)\n";
    std::cout << "  \033[1;36minterpolate\033[0m   Interpolate dataset (finite differences, Lagrange, splines)\n";
    std::cout << "  \033[1;36msolve-system\033[0m  Solve linear system Ax = b with row-operation tracking\n";
    std::cout << "  \033[1;36mfactorize\033[0m     Matrix factorizations (Doolittle, Crout, Cholesky)\n";
    std::cout << "  \033[1;36meigen\033[0m         Dominant eigenvalue & eigenvector (Power Method)\n";
    std::cout << "  \033[1;36mspline\033[0m        Cubic spline & B-spline curve approximation\n";
    std::cout << "  \033[1;36mode\033[0m           Solve initial value ODE dy/dx = f(x, y)\n";
    std::cout << "  \033[1;36mpde\033[0m           Solve partial differential equations (Heat, Laplace, Wave)\n";
    std::cout << "  \033[1;36mtaylor\033[0m        Taylor series polynomial expansion & error analysis\n";
    std::cout << "  \033[1;36mcompare\033[0m       Multi-method benchmark & error comparison\n";
    std::cout << "  \033[1;36madvise\033[0m        Heuristic Numerical Method Advisor\n";
    std::cout << "  \033[1;36mproblem\033[0m       Curated engineering problems (Mechanical, Civil, Electrical, Thermal)\n";
    std::cout << "  \033[1;36marithmetic\033[0m    Floating-point stability & catastrophic cancellation lab\n\n";

    std::cout << "\033[1mEXAMPLES:\033[0m\n";
    std::cout << "  ees root \"x^3 - 4*x - 9\" --method newton --x0 2.5 --plot\n";
    std::cout << "  ees root \"x^3 - 4*x - 9\" --compare\n";
    std::cout << "  ees optimize \"(x - 3)^2 + 4\" --a 0 --b 5 --method golden\n";
    std::cout << "  ees integrate \"4 / (1 + x^2)\" --from 0 --to 1 --method gauss5\n";
    std::cout << "  ees solve-system --matrix \"[[2,1,-1],[-3,-1,2],[-2,1,2]]\" --rhs \"[8,-11,-3]\" --steps\n";
    std::cout << "  ees factorize --matrix \"[[4,12],[12,45]]\" --method cholesky\n";
    std::cout << "  ees eigen --matrix \"[[2,2],[1,3]]\"\n";
    std::cout << "  ees ode \"y\" --x0 0 --y0 1 --xEnd 1 --h 0.1 --method rk4\n";
    std::cout << "  ees pde heat --alpha 0.1 --L 1.0 --tMax 0.5\n";
    std::cout << "  ees demo\n";
}

bool CLIApp::hasFlag(const std::vector<std::string>& args, const std::string& flag) {
    for (const auto& a : args) {
        if (a == flag) return true;
    }
    return false;
}

std::string CLIApp::getOption(const std::vector<std::string>& args, const std::string& opt, const std::string& defaultVal) {
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == opt && i + 1 < args.size()) {
            return args[i + 1];
        }
    }
    return defaultVal;
}

std::vector<double> CLIApp::parseCSV(const std::string& str) {
    std::vector<double> res;
    std::string clean = str;
    clean.erase(std::remove(clean.begin(), clean.end(), '['), clean.end());
    clean.erase(std::remove(clean.begin(), clean.end(), ']'), clean.end());

    std::stringstream ss(clean);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) {
            res.push_back(std::stod(item));
        }
    }
    return res;
}

std::vector<std::vector<double>> CLIApp::parseMatrix(const std::string& str) {
    std::vector<std::vector<double>> mat;
    size_t i = 0;
    while (i < str.length()) {
        if (str[i] == '[') {
            size_t end = str.find(']', i + 1);
            if (end != std::string::npos) {
                std::string rowStr = str.substr(i + 1, end - i - 1);
                auto rowVals = parseCSV(rowStr);
                if (!rowVals.empty()) {
                    mat.push_back(rowVals);
                }
                i = end + 1;
                continue;
            }
        }
        ++i;
    }
    return mat;
}

int CLIApp::handleRoot(const std::vector<std::string>& args) {
    if (args.empty()) {
        std::cerr << "Error: Expression required for root-finding.\n";
        return 1;
    }
    std::string expr = args[0];
    core::ExpressionParser parser(expr);
    if (!parser.isValid()) {
        std::cerr << "Invalid expression: " << parser.getErrorMessage() << "\n";
        return 1;
    }
    auto f = [parser](double x) mutable { return parser.evaluate(x); };

    double a = std::stod(getOption(args, "--a", "0.0"));
    double b = std::stod(getOption(args, "--b", "3.0"));
    double x0 = std::stod(getOption(args, "--x0", std::to_string(0.5 * (a + b))));
    double tol = std::stod(getOption(args, "--tol", "1e-7"));
    int maxIter = std::stoi(getOption(args, "--max-iter", "100"));

    if (hasFlag(args, "--compare")) {
        std::cout << unit1::RootFindingLab::compareRootMethods(f, a, b, x0, tol, maxIter) << "\n";
        return 0;
    }

    std::string method = getOption(args, "--method", "newton");
    core::NumericalResult res;
    if (method == "bisection") {
        res = unit1::RootFindingLab::bisection(f, a, b, tol, maxIter);
    } else if (method == "secant") {
        res = unit1::RootFindingLab::secant(f, a, b, tol, maxIter);
    } else {
        res = unit1::RootFindingLab::newtonRaphson(f, x0, tol, maxIter);
    }

    std::cout << res.summary() << "\n";
    if (hasFlag(args, "--step") || hasFlag(args, "--table")) {
        std::cout << res.tableString({"f'(x)", "a", "b", "interval_width"}) << "\n";
    }
    if (hasFlag(args, "--plot")) {
        std::cout << Visualizer::plotFunction(f, a - 0.5, b + 0.5, 65, 15, {res.solution}, "Root Finding Plot") << "\n";
        std::cout << Visualizer::plotConvergence(res.history, 65, 14) << "\n";
    }
    return res.success ? 0 : 1;
}

int CLIApp::handleOptimize(const std::vector<std::string>& args) {
    if (args.empty()) {
        std::cerr << "Error: Expression required for optimization.\n";
        return 1;
    }
    std::string expr = args[0];
    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };

    double a = std::stod(getOption(args, "--a", "0.0"));
    double b = std::stod(getOption(args, "--b", "5.0"));
    double x0 = std::stod(getOption(args, "--x0", "1.0"));
    double tol = std::stod(getOption(args, "--tol", "1e-5"));

    std::string method = getOption(args, "--method", "golden");
    core::NumericalResult res;
    if (method == "fibonacci") {
        res = unit1::OptimizationLab::fibonacciSearch(f, a, b, tol);
    } else if (method == "newton") {
        res = unit1::OptimizationLab::newtonsMethod1D(f, x0, tol);
    } else {
        res = unit1::OptimizationLab::goldenSectionSearch(f, a, b, tol);
    }

    std::cout << res.summary() << "\n";
    if (hasFlag(args, "--step") || hasFlag(args, "--table")) {
        std::cout << res.tableString({"a", "b", "width"}) << "\n";
    }
    if (hasFlag(args, "--plot")) {
        std::cout << Visualizer::plotIntervalShrinkage(res.history, 50) << "\n";
    }
    return 0;
}

int CLIApp::handleIntegrate(const std::vector<std::string>& args) {
    if (args.empty()) {
        std::cerr << "Error: Expression required for integration.\n";
        return 1;
    }
    std::string expr = args[0];
    core::ExpressionParser parser(expr);
    auto f = [parser](double x) mutable { return parser.evaluate(x); };

    double a = std::stod(getOption(args, "--from", "0.0"));
    double b = std::stod(getOption(args, "--to", "1.0"));
    int n = std::stoi(getOption(args, "--n", "100"));
    double trueVal = std::stod(getOption(args, "--exact", "0.0"));

    if (hasFlag(args, "--compare")) {
        std::cout << unit2::IntegrationLab::compareIntegrationMethods(f, a, b, trueVal) << "\n";
        return 0;
    }

    std::string method = getOption(args, "--method", "simpson13");
    unit2::IntegrationResult res;
    if (method == "trapezoidal") {
        res = unit2::IntegrationLab::trapezoidal(f, a, b, n, trueVal);
    } else if (method == "simpson38") {
        res = unit2::IntegrationLab::simpson38(f, a, b, n, trueVal);
    } else if (method == "romberg") {
        res = unit2::IntegrationLab::romberg(f, a, b, 5, 1e-9, trueVal);
    } else if (method == "gauss2") {
        res = unit2::IntegrationLab::gaussianQuadrature(f, a, b, 2, trueVal);
    } else if (method == "gauss3") {
        res = unit2::IntegrationLab::gaussianQuadrature(f, a, b, 3, trueVal);
    } else if (method == "gauss4") {
        res = unit2::IntegrationLab::gaussianQuadrature(f, a, b, 4, trueVal);
    } else if (method == "gauss5") {
        res = unit2::IntegrationLab::gaussianQuadrature(f, a, b, 5, trueVal);
    } else {
        res = unit2::IntegrationLab::simpson13(f, a, b, n, trueVal);
    }

    std::cout << res.summary() << "\n";
    return 0;
}

int CLIApp::handleInterpolate(const std::vector<std::string>& args) {
    std::string xStr = getOption(args, "--x", "1,2,3,4");
    std::string yStr = getOption(args, "--y", "0,5,22,57");
    double evalX = std::stod(getOption(args, "--eval", "2.5"));

    core::Vector x = parseCSV(xStr);
    core::Vector y = parseCSV(yStr);

    if (x.size() != y.size() || x.empty()) {
        std::cerr << "Error: Equal number of x and y points required.\n";
        return 1;
    }

    if (hasFlag(args, "--diff-table")) {
        auto dt = unit2::InterpolationLab::computeForwardDifferenceTable(x, y);
        std::cout << dt.toString() << "\n";
    }

    if (hasFlag(args, "--compare")) {
        std::cout << unit2::InterpolationLab::compareInterpolation(x, y, evalX) << "\n";
        return 0;
    }

    std::string method = getOption(args, "--method", "newton-forward");
    unit2::InterpolationResult res;
    if (method == "backward" || method == "newton-backward") {
        res = unit2::InterpolationLab::gregoryNewtonBackward(x, y, evalX);
    } else if (method == "lagrange") {
        res = unit2::InterpolationLab::lagrange(x, y, evalX);
    } else if (method == "divided") {
        res = unit2::InterpolationLab::newtonDividedDifference(x, y, evalX);
    } else {
        res = unit2::InterpolationLab::gregoryNewtonForward(x, y, evalX);
    }

    std::cout << res.summary() << "\n";
    return 0;
}

int CLIApp::handleLinearSystem(const std::vector<std::string>& args) {
    std::string matStr = getOption(args, "--matrix", "[[2,1,-1],[-3,-1,2],[-2,1,2]]");
    std::string rhsStr = getOption(args, "--rhs", "[8,-11,-3]");

    auto matData = parseMatrix(matStr);
    auto bVec = parseCSV(rhsStr);

    if (matData.empty() || matData.size() != bVec.size()) {
        std::cerr << "Error: Incompatible matrix and RHS dimensions.\n";
        return 1;
    }

    core::Matrix A(matData.size(), matData[0].size());
    for (size_t r = 0; r < matData.size(); ++r) {
        for (size_t c = 0; c < matData[0].size(); ++c) {
            A(r, c) = matData[r][c];
        }
    }
    core::Vector b = bVec;

    std::string method = getOption(args, "--method", "gauss");
    unit3::LinearSystemResult res;
    if (method == "jordan") {
        res = unit3::LinearAlgebraLab::gaussJordan(A, b);
    } else {
        res = unit3::LinearAlgebraLab::gaussElimination(A, b, true);
    }

    std::cout << res.summary() << "\n";
    if (hasFlag(args, "--steps")) {
        std::cout << "\nElimination Step-by-Step Operations:\n";
        for (const auto& step : res.stepHistory) {
            std::cout << "  [" << step.stepNumber << "] " << step.operation << "\n";
        }
    }
    return res.success ? 0 : 1;
}

int CLIApp::handleFactorize(const std::vector<std::string>& args) {
    std::string matStr = getOption(args, "--matrix", "[[4,12],[12,45]]");
    auto matData = parseMatrix(matStr);
    if (matData.empty()) {
        std::cerr << "Error: Valid square matrix required.\n";
        return 1;
    }

    core::Matrix A(matData.size(), matData[0].size());
    for (size_t r = 0; r < matData.size(); ++r) {
        for (size_t c = 0; c < matData[0].size(); ++c) {
            A(r, c) = matData[r][c];
        }
    }

    std::string method = getOption(args, "--method", "cholesky");
    unit3::FactorizationResult res;
    if (method == "doolittle") {
        res = unit3::LinearAlgebraLab::doolittleLU(A);
    } else if (method == "crout") {
        res = unit3::LinearAlgebraLab::croutLU(A);
    } else {
        res = unit3::LinearAlgebraLab::cholesky(A);
    }

    std::cout << res.summary() << "\n";
    return res.success ? 0 : 1;
}

int CLIApp::handleEigen(const std::vector<std::string>& args) {
    std::string matStr = getOption(args, "--matrix", "[[2,2],[1,3]]");
    auto matData = parseMatrix(matStr);
    if (matData.empty()) {
        std::cerr << "Error: Valid square matrix required.\n";
        return 1;
    }

    core::Matrix A(matData.size(), matData[0].size());
    for (size_t r = 0; r < matData.size(); ++r) {
        for (size_t c = 0; c < matData[0].size(); ++c) {
            A(r, c) = matData[r][c];
        }
    }

    double tol = std::stod(getOption(args, "--tol", "1e-7"));
    auto res = unit3::LinearAlgebraLab::powerMethod(A, core::Vector(), tol);
    std::cout << res.summary() << "\n";
    if (hasFlag(args, "--plot")) {
        std::cout << Visualizer::plotConvergence(res.history, 65, 14, "Power Method Convergence") << "\n";
    }
    return res.success ? 0 : 1;
}

int CLIApp::handleSpline(const std::vector<std::string>& args) {
    std::string xStr = getOption(args, "--x", "0,1,2,3");
    std::string yStr = getOption(args, "--y", "0,1,8,27");
    double evalX = std::stod(getOption(args, "--eval", "1.5"));

    core::Vector x = parseCSV(xStr);
    core::Vector y = parseCSV(yStr);

    std::string method = getOption(args, "--method", "cubic");
    unit3::SplineResult res;
    if (method == "linear") {
        res = unit3::SplineLab::linearSpline(x, y, evalX);
    } else if (method == "quadratic") {
        res = unit3::SplineLab::quadraticSpline(x, y, evalX);
    } else {
        res = unit3::SplineLab::naturalCubicSpline(x, y, evalX);
    }

    std::cout << res.summary() << "\n";
    return 0;
}

int CLIApp::handleODE(const std::vector<std::string>& args) {
    if (args.empty()) {
        std::cerr << "Error: Expression required for ODE dy/dx = f(x, y).\n";
        return 1;
    }
    std::string expr = args[0];
    core::ExpressionParser parser(expr);
    auto f = [parser](double x, double y) mutable {
        return parser.evaluate(x, y);
    };

    double x0 = std::stod(getOption(args, "--x0", "0.0"));
    double y0 = std::stod(getOption(args, "--y0", "1.0"));
    double xEnd = std::stod(getOption(args, "--xEnd", "1.0"));
    double h = std::stod(getOption(args, "--h", "0.1"));

    if (hasFlag(args, "--compare")) {
        std::cout << unit4::ODESolverLab::compareODESolvers(f, x0, y0, xEnd, h) << "\n";
        return 0;
    }

    std::string method = getOption(args, "--method", "rk4");
    unit4::ODESolution res;
    if (method == "euler") {
        res = unit4::ODESolverLab::euler(f, x0, y0, xEnd, h);
    } else if (method == "heun") {
        res = unit4::ODESolverLab::eulerPredictorCorrector(f, x0, y0, xEnd, h);
    } else if (method == "adams") {
        res = unit4::ODESolverLab::adamsBashforth(f, x0, y0, xEnd, h);
    } else if (method == "milne") {
        res = unit4::ODESolverLab::milne(f, x0, y0, xEnd, h);
    } else {
        res = unit4::ODESolverLab::rungeKutta4(f, x0, y0, xEnd, h);
    }

    std::cout << res.summary() << "\n";
    if (hasFlag(args, "--table")) {
        std::cout << res.tableString(15) << "\n";
    }
    return 0;
}

int CLIApp::handlePDE(const std::vector<std::string>& args) {
    if (args.empty()) {
        std::cerr << "Error: PDE type required: 'heat', 'laplace', or 'wave'.\n";
        return 1;
    }
    std::string type = args[0];
    unit4::PDESolution sol;

    if (type == "laplace") {
        sol = unit4::PDESolverLab::solveLaplaceEquation2D(1.0, 1.0, 15, 15, 100.0, 0.0, 0.0, 0.0, 1.3, 1e-4, 500);
    } else if (type == "wave") {
        auto u0 = [](double x){ return std::sin(3.141592653589793 * x); };
        sol = unit4::PDESolverLab::solveWaveEquation1D(1.0, 1.0, 1.0, 21, 51, u0, nullptr, 0.0, 0.0);
    } else {
        double alpha = std::stod(getOption(args, "--alpha", "0.05"));
        double L = std::stod(getOption(args, "--L", "1.0"));
        double tMax = std::stod(getOption(args, "--tMax", "0.5"));
        auto u0 = [](double x){ return 100.0 * std::sin(3.141592653589793 * x); };
        sol = unit4::PDESolverLab::solveHeatEquation1D(alpha, L, tMax, 21, 51, u0, 0.0, 0.0);
    }

    std::cout << sol.summary() << "\n";
    std::cout << sol.asciiHeatmap(45, 12) << "\n";
    return 0;
}

int CLIApp::handleTaylor(const std::vector<std::string>& args) {
    std::string expr = args.empty() ? "exp(x)" : args[0];
    double x0 = std::stod(getOption(args, "--x0", "0.0"));
    double x = std::stod(getOption(args, "--x", "1.0"));
    int order = std::stoi(getOption(args, "--order", "5"));
    double trueVal = std::stod(getOption(args, "--exact", std::to_string(std::exp(x))));

    auto res = unit1::FoundationsLab::computeTaylorSeriesExpr(expr, x0, x, order, trueVal);
    std::cout << res.summary() << "\n";
    return 0;
}

int CLIApp::handleAdvise(const std::vector<std::string>& args) {
    advisor::ProblemProfile prof;
    prof.category = advisor::ProblemCategory::RootFinding;

    if (!args.empty()) {
        std::string cat = args[0];
        if (cat == "linear") prof.category = advisor::ProblemCategory::LinearSystem;
        else if (cat == "integration") prof.category = advisor::ProblemCategory::NumericalIntegration;
        else if (cat == "interpolation") prof.category = advisor::ProblemCategory::Interpolation;
        else if (cat == "ode") prof.category = advisor::ProblemCategory::ODEInitialValue;
    }

    prof.hasBracketingInterval = hasFlag(args, "--bracket");
    prof.hasAnalyticDerivative = hasFlag(args, "--derivative");
    prof.isMatrixSymmetricPositiveDefinite = hasFlag(args, "--spd");

    auto rec = advisor::MethodAdvisor::advise(prof);
    std::cout << rec.toString() << "\n";
    return 0;
}

int CLIApp::handleProblem(const std::vector<std::string>& args) {
    if (args.empty()) {
        auto problems = problems::ProblemLibrary::getAllProblems();
        std::cout << "Curated Engineering Problems:\n";
        for (const auto& p : problems) {
            std::cout << "  " << p.id << " : " << p.title << "\n";
        }
        return 0;
    }

    std::string id = args[0];
    try {
        auto p = problems::ProblemLibrary::getProblemById(id);
        std::cout << "Problem: " << p.title << " (" << p.id << ")\n";
        std::cout << "Model:   " << p.mathematicalModel << "\n";
        std::cout << "Result:  " << p.expectedResult << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}

int CLIApp::execute(int argc, char* argv[]) {
    if (argc < 2) {
        // If no arguments provided, launch the interactive Workbench directly!
        Workbench workbench;
        workbench.run();
        return 0;
    }

    std::string cmd = argv[1];

    if (cmd == "--help" || cmd == "-h" || cmd == "help") {
        printUsage();
        return 0;
    }
    if (cmd == "--version" || cmd == "-v" || cmd == "version") {
        printVersion();
        return 0;
    }
    if (cmd == "demo" || cmd == "--demo" || cmd == "-d") {
        Workbench workbench;
        workbench.runDemoMode(false);
        return 0;
    }
    if (cmd == "interactive" || cmd == "--interactive" || cmd == "-i") {
        Workbench workbench;
        workbench.run();
        return 0;
    }
    if (cmd == "arithmetic") {
        std::cout << unit1::FoundationsLab::runArithmeticStabilityLab() << "\n";
        return 0;
    }

    std::vector<std::string> subArgs;
    for (int i = 2; i < argc; ++i) {
        subArgs.push_back(argv[i]);
    }

    if (cmd == "root")         return handleRoot(subArgs);
    if (cmd == "optimize")     return handleOptimize(subArgs);
    if (cmd == "integrate")    return handleIntegrate(subArgs);
    if (cmd == "interpolate")  return handleInterpolate(subArgs);
    if (cmd == "solve-system") return handleLinearSystem(subArgs);
    if (cmd == "factorize")    return handleFactorize(subArgs);
    if (cmd == "eigen")        return handleEigen(subArgs);
    if (cmd == "spline")       return handleSpline(subArgs);
    if (cmd == "ode")          return handleODE(subArgs);
    if (cmd == "pde")          return handlePDE(subArgs);
    if (cmd == "taylor")       return handleTaylor(subArgs);
    if (cmd == "advise")       return handleAdvise(subArgs);
    if (cmd == "problem")      return handleProblem(subArgs);

    std::cerr << "Unknown command: '" << cmd << "'. Run 'ees --help' for available commands.\n";
    return 1;
}

} // namespace ees::ui

