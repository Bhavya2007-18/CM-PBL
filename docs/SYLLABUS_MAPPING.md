# SYLLABUS MAPPING — COMPUTATIONAL METHODS PBL

This document proves exhaustive 100% curriculum coverage mapping every syllabus topic directly to the manual C++ implementation files, header abstractions, automated test cases, and workbench demonstration options.

| Syllabus Topic | Application Module | Header File | Implementation File | Test & Verification | Workbench Demonstration |
|---|---|---|---|---|---|
| **UNIT I: Numerical Foundations** |
| Taylor Series | `ees::unit1::FoundationsLab` | `include/ees/unit1/Foundations.hpp` | `src/unit1/Foundations.cpp` | `tests/test_unit1.cpp` | Dashboard Option 1 -> Sub-option 2 |
| Rolle's Theorem | `ees::unit1::FoundationsLab` | `include/ees/unit1/Foundations.hpp` | `src/unit1/Foundations.cpp` | `tests/test_unit1.cpp` | Dashboard Option 1 -> Sub-option 3 |
| Mean Value Theorem | `ees::unit1::FoundationsLab` | `include/ees/unit1/Foundations.hpp` | `src/unit1/Foundations.cpp` | `tests/test_unit1.cpp` | Dashboard Option 1 -> Sub-option 3 |
| Approximations & Error Analysis | `ees::core::ErrorAnalysis` | `include/ees/core/ErrorAnalysis.hpp` | `src/core/ErrorAnalysis.cpp` | `tests/test_core.cpp` | Dashboard Option 1 -> Sub-option 1 |
| Absolute / Relative / % Error | `ees::core::ErrorAnalysis` | `include/ees/core/ErrorAnalysis.hpp` | `src/core/ErrorAnalysis.cpp` | `tests/test_core.cpp` | Formatted in all method tables |
| Significant Digits (Scarborough) | `ees::core::ErrorAnalysis` | `include/ees/core/ErrorAnalysis.hpp` | `src/core/ErrorAnalysis.cpp` | `tests/test_core.cpp` | Dashboard Option 1 -> Sub-option 1 |
| Data Representation & IEEE-754 | `ees::core::ErrorAnalysis` | `include/ees/core/ErrorAnalysis.hpp` | `src/core/ErrorAnalysis.cpp` | `tests/test_core.cpp` | Dashboard Option 1 -> Sub-option 1 |
| Computer Arithmetic & Round-off | `ees::core::ErrorAnalysis` | `include/ees/core/ErrorAnalysis.hpp` | `src/core/ErrorAnalysis.cpp` | `tests/test_core.cpp` | Dashboard Option 1 -> Sub-option 1 |
| Loss of Significance / Cancellation | `ees::unit1::FoundationsLab` | `include/ees/unit1/Foundations.hpp` | `src/unit1/Foundations.cpp` | `tests/test_unit1.cpp` | Dashboard Option 1 -> Sub-option 1 |
| **UNIT I: Root Location** |
| Bisection Method | `ees::unit1::RootFindingLab` | `include/ees/unit1/RootFinding.hpp` | `src/unit1/RootFinding.cpp` | `tests/test_unit1.cpp` | Dashboard Option 2 & Option 14 [1/9] |
| Newton-Raphson Method | `ees::unit1::RootFindingLab` | `include/ees/unit1/RootFinding.hpp` | `src/unit1/RootFinding.cpp` | `tests/test_unit1.cpp` | Dashboard Option 2 & Option 14 [1/9] |
| Secant Method | `ees::unit1::RootFindingLab` | `include/ees/unit1/RootFinding.hpp` | `src/unit1/RootFinding.cpp` | `tests/test_unit1.cpp` | Dashboard Option 2 |
| Convergence Rate Analysis | `ees::core::NumericalResult` | `include/ees/core/NumericalResult.hpp` | `src/core/NumericalResult.cpp` | `tests/test_core.cpp` | Convergence rate classification in all solvers |
| Iteration Tracking & Trajectories | `ees::core::IterationRecord` | `include/ees/core/NumericalResult.hpp` | `src/core/NumericalResult.cpp` | `tests/test_core.cpp` | Table String output across all modules |
| **UNIT I: Optimization** |
| Fibonacci Search | `ees::unit1::OptimizationLab` | `include/ees/unit1/Optimization.hpp` | `src/unit1/Optimization.cpp` | `tests/test_unit1.cpp` | Dashboard Option 3 -> Sub-option 1 |
| Golden Section Search | `ees::unit1::OptimizationLab` | `include/ees/unit1/Optimization.hpp` | `src/unit1/Optimization.cpp` | `tests/test_unit1.cpp` | Dashboard Option 3 -> Sub-option 1 |
| Newton's Method (1D Opt) | `ees::unit1::OptimizationLab` | `include/ees/unit1/Optimization.hpp` | `src/unit1/Optimization.cpp` | `tests/test_unit1.cpp` | Dashboard Option 3 -> Sub-option 1 |
| Steepest Descent (Multivariate) | `ees::unit1::OptimizationLab` | `include/ees/unit1/Optimization.hpp` | `src/unit1/Optimization.cpp` | `tests/test_unit1.cpp` | Dashboard Option 3 -> Sub-option 2 |
| Nelder-Mead Simplex | `ees::unit1::OptimizationLab` | `include/ees/unit1/Optimization.hpp` | `src/unit1/Optimization.cpp` | `tests/test_unit1.cpp` | Dashboard Option 3 & Option 14 [2/9] |
| **UNIT II: Interpolation** |
| Finite Difference Tables | `ees::unit2::InterpolationLab` | `include/ees/unit2/Interpolation.hpp` | `src/unit2/Interpolation.cpp` | `tests/test_unit2.cpp` | Dashboard Option 4 |
| Gregory-Newton Forward | `ees::unit2::InterpolationLab` | `include/ees/unit2/Interpolation.hpp` | `src/unit2/Interpolation.cpp` | `tests/test_unit2.cpp` | Dashboard Option 4 & Option 14 [3/9] |
| Gregory-Newton Backward | `ees::unit2::InterpolationLab` | `include/ees/unit2/Interpolation.hpp` | `src/unit2/Interpolation.cpp` | `tests/test_unit2.cpp` | Dashboard Option 4 |
| Lagrange Interpolation | `ees::unit2::InterpolationLab` | `include/ees/unit2/Interpolation.hpp` | `src/unit2/Interpolation.cpp` | `tests/test_unit2.cpp` | Dashboard Option 4 |
| Newton Divided Difference | `ees::unit2::InterpolationLab` | `include/ees/unit2/Interpolation.hpp` | `src/unit2/Interpolation.cpp` | `tests/test_unit2.cpp` | Dashboard Option 4 |
| Errors in Polynomial Interpolation | `ees::unit2::InterpolationLab` | `include/ees/unit2/Interpolation.hpp` | `src/unit2/Interpolation.cpp` | `tests/test_unit2.cpp` | Displayed in polynomial summary |
| **UNIT II: Numerical Integration** |
| Definite Integrals | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 |
| Newton-Cotes Quadrature | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 |
| Trapezoidal Rule | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 & Option 11 |
| Simpson's 1/3 Rule | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 & Option 11 |
| Simpson's 3/8 Rule | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 & Option 11 |
| Romberg Algorithm | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 |
| Gaussian Quadrature (Legendre) | `ees::unit2::IntegrationLab` | `include/ees/unit2/Integration.hpp` | `src/unit2/Integration.cpp` | `tests/test_unit2.cpp` | Dashboard Option 5 & Option 14 [4/9] |
| **UNIT III: Linear Systems** |
| Existence & Determinant | `ees::core::Matrix` | `include/ees/core/Matrix.hpp` | `src/core/Matrix.cpp` | `tests/test_core.cpp` | Dashboard Option 6 |
| Gauss Elimination | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 & Option 14 [5/9] |
| Partial Pivoting Lab | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 -> Sub-option 4 |
| Gauss-Jordan & RREF | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 -> Sub-option 2 |
| Matrix Inversion | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 -> Sub-option 2 |
| Doolittle LU Factorization | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 -> Sub-option 3 |
| Crout LU Factorization | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 -> Sub-option 3 |
| Cholesky LL^T Factorization | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 6 -> Sub-option 3 |
| Power Method (Eigenvalues) | `ees::unit3::LinearAlgebraLab` | `include/ees/unit3/LinearSystems.hpp` | `src/unit3/LinearSystems.cpp` | `tests/test_unit3.cpp` | Dashboard Option 7 & Option 14 [6/9] |
| First-Degree Spline | `ees::unit3::SplineLab` | `include/ees/unit3/Splines.hpp` | `src/unit3/Splines.cpp` | `tests/test_unit3.cpp` | Dashboard Option 8 |
| Second-Degree Spline | `ees::unit3::SplineLab` | `include/ees/unit3/Splines.hpp` | `src/unit3/Splines.cpp` | `tests/test_unit3.cpp` | Dashboard Option 8 |
| Natural Cubic Splines | `ees::unit3::SplineLab` | `include/ees/unit3/Splines.hpp` | `src/unit3/Splines.cpp` | `tests/test_unit3.cpp` | Dashboard Option 8 & Option 14 [7/9] |
| B-Splines (Cox-de Boor) | `ees::unit3::SplineLab` | `include/ees/unit3/Splines.hpp` | `src/unit3/Splines.cpp` | `tests/test_unit3.cpp` | Dashboard Option 8 |
| Interpolation vs Approximation | `ees::unit3::SplineLab` | `include/ees/unit3/Splines.hpp` | `src/unit3/Splines.cpp` | `tests/test_unit3.cpp` | Displayed in Spline summary |
| **UNIT IV: Differential Equations** |
| Picard's Method | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | Verified in unit4 test | Dashboard Option 9 |
| Taylor Series Method ODE | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | Verified in unit4 test | Dashboard Option 9 |
| Euler Method | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 9 & Option 11 |
| Runge-Kutta 4th-Order (RK4) | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 9 & Option 14 [8/9] |
| Euler Predictor-Corrector (Heun) | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 9 & Option 11 |
| Adams-Bashforth-Moulton | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 9 & Option 11 |
| Milne's Predictor-Corrector | `ees::unit4::ODESolverLab` | `include/ees/unit4/ODESolvers.hpp` | `src/unit4/ODESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 9 & Option 11 |
| Parabolic PDE (Heat Equation) | `ees::unit4::PDESolverLab` | `include/ees/unit4/PDESolvers.hpp` | `src/unit4/PDESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 10 & Option 14 [9/9] |
| Elliptic PDE (Laplace / SOR) | `ees::unit4::PDESolverLab` | `include/ees/unit4/PDESolvers.hpp` | `src/unit4/PDESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 10 |
| Hyperbolic PDE (Wave Equation) | `ees::unit4::PDESolverLab` | `include/ees/unit4/PDESolvers.hpp` | `src/unit4/PDESolvers.cpp` | `tests/test_unit4.cpp` | Dashboard Option 10 |
| **Workbench Features** |
| Method Comparison Engine | `ees::comparison::MethodComparator` | `include/ees/comparison/MethodComparator.hpp` | `src/comparison/MethodComparator.cpp` | Integrated across all units | Dashboard Option 11 |
| Engineering Problem Library | `ees::problems::ProblemLibrary` | `include/ees/problems/ProblemLibrary.hpp` | `src/problems/ProblemLibrary.cpp` | Verified with exact answers | Dashboard Option 12 |
| Numerical Method Advisor | `ees::advisor::MethodAdvisor` | `include/ees/advisor/MethodAdvisor.hpp` | `src/advisor/MethodAdvisor.cpp` | Verified heuristic rules | Dashboard Option 13 |
| One-Click PBL Viva Demo Mode | `ees::ui::Workbench` | `include/ees/ui/Workbench.hpp` | `src/ui/Workbench.cpp` | Verified via `--demo` | Dashboard Option 14 |

