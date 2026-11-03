# EES — Engineering Equation Solver
### Numerical Methods & Engineering Analysis Workbench
**A Unified Computational Methods / Numerical Methods PBL Platform (C++17)**

![C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Build-Tool](https://img.shields.io/badge/Build-CMake%20%7C%20Ninja-success.svg)
![Tests](https://img.shields.io/badge/Tests-16%2F16%20Passing%20(100%25)-brightgreen.svg)
![Curriculum](https://img.shields.io/badge/Curriculum-Units%20I--IV%20(100%25)-orange.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%2010%2F11-blueviolet.svg)

---

## 1. Project Identity & Philosophy

> **Solve → Explain → Visualize → Compare → Understand**

**EES (Engineering Equation Solver)** is an engineering-grade scientific computing application engineered from first principles in modern C++ (C++17/20). Rather than a disconnected set of numerical scripts, EES provides an integrated **Numerical Methods Engineering Workbench** designed for university Computational Methods PBL (Project-Based Learning) demonstrations, viva defense, and engineering design.

### Key Highlights
- **100% Pure C++ Implementation**: All algorithms (roots, optimization, matrices, eigenvalues, splines, ODEs, PDEs) are implemented from scratch without delegating to Python, MATLAB, NumPy, SciPy, or external numerical black boxes.
- **Layered Architecture**: Complete isolation between Core Math (`Vector`, `Matrix`, `ExpressionParser`, `ErrorAnalysis`), Numerical Engines (Units I–IV), Application Controller (Comparison, Advisor, Problems), and Visualization/UI.
- **Scientific Visualization Engine**: Generates 2D function curves, semi-log convergence trajectories, interval shrinkage bars, and PDE spatial-temporal heatmaps directly in the Windows terminal, plus standalone HTML5 dashboard export.
- **Step-by-Step & Method Comparison**: Inspect every iteration, row operation, and convergence step; benchmark competing methods on identical problems.
- **Direct Academic Syllabus Mapping**: Directly maps Units I, II, III, and IV of standard B.Tech Computational Methods curriculum.

---

## 2. Academic Syllabus Mapping Matrix

| Unit | Syllabus Topics | Implemented Algorithms | C++ Implementation Files |
|---|---|---|---|
| **UNIT I** | Numerical Foundations | Taylor Series, Rolle's Theorem, Mean Value Theorem, Absolute/Relative/% Error, Scarborough Significant Digits, IEEE-754 Bit Inspection, Catastrophic Cancellation | [`src/core/ErrorAnalysis.cpp`](file:///d:/CODING/CM%20PBL/src/core/ErrorAnalysis.cpp)<br>[`src/unit1/Foundations.cpp`](file:///d:/CODING/CM%20PBL/src/unit1/Foundations.cpp) |
| **UNIT I** | Root Location | Bisection Method (Linear, Order 1), Newton-Raphson Method (Quadratic, Order 2), Secant Method (Superlinear, Order 1.62), Iteration Tracking | [`src/unit1/RootFinding.cpp`](file:///d:/CODING/CM%20PBL/src/unit1/RootFinding.cpp) |
| **UNIT I** | Optimization | Fibonacci Search, Golden Section Search (0.618 ratio), Newton's 1D Optimization, Steepest Descent with Armijo Backtracking, Nelder-Mead Simplex Algorithm | [`src/unit1/Optimization.cpp`](file:///d:/CODING/CM%20PBL/src/unit1/Optimization.cpp) |
| **UNIT II** | Interpolation | Forward/Backward/Divided Difference Tables, Gregory-Newton Forward, Gregory-Newton Backward, Lagrange Polynomial, Newton Divided Difference | [`src/unit2/Interpolation.cpp`](file:///d:/CODING/CM%20PBL/src/unit2/Interpolation.cpp) |
| **UNIT II** | Numerical Integration | Composite Trapezoidal ($O(h^2)$), Simpson's 1/3 ($O(h^4)$), Simpson's 3/8 ($O(h^4)$), Romberg Extrapolation Table, Gauss-Legendre Quadrature (2, 3, 4, 5-point exact nodes/weights) | [`src/unit2/Integration.cpp`](file:///d:/CODING/CM%20PBL/src/unit2/Integration.cpp) |
| **UNIT III** | Linear Systems | Gauss Elimination with Partial Pivoting ($2n^3/3$ FLOPs), Pivoting Stability Lab, Gauss-Jordan (Reduced Row Echelon Form, Matrix Inversion) | [`src/unit3/LinearSystems.cpp`](file:///d:/CODING/CM%20PBL/src/unit3/LinearSystems.cpp) |
| **UNIT III** | Matrix Factorization | Doolittle LU ($L_{ii}=1$), Crout LU ($U_{ii}=1$), Cholesky $A = LL^T$ (Symmetric Positive Definite, $n^3/3$ FLOPs, Frobenius reconstruction verification) | [`src/unit3/LinearSystems.cpp`](file:///d:/CODING/CM%20PBL/src/unit3/LinearSystems.cpp) |
| **UNIT III** | Eigenvalues & Splines | Power Method with Rayleigh Quotient Iteration, Linear Spline, Quadratic Spline, Natural Cubic Spline ($C^2$ continuous, Thomas Tridiagonal Solver), B-Spline (Cox-de Boor) | [`src/unit3/LinearSystems.cpp`](file:///d:/CODING/CM%20PBL/src/unit3/LinearSystems.cpp)<br>[`src/unit3/Splines.cpp`](file:///d:/CODING/CM%20PBL/src/unit3/Splines.cpp) |
| **UNIT IV** | Numerical ODEs | Picard's Method, Taylor Series Method, Euler Method ($O(h)$), Heun Predictor-Corrector ($O(h^2)$), Classic Runge-Kutta 4th-Order (RK4, $O(h^4)$), Adams-Bashforth-Moulton 4-step, Milne's Predictor-Corrector | [`src/unit4/ODESolvers.cpp`](file:///d:/CODING/CM%20PBL/src/unit4/ODESolvers.cpp) |
| **UNIT IV** | Numerical PDEs | Parabolic 1D Heat Equation (FTCS scheme with CFL condition $r \le 0.5$), Elliptic 2D Laplace Equation (Gauss-Seidel with SOR acceleration), Hyperbolic 1D Wave Equation (Central differences, Courant $C \le 1.0$) | [`src/unit4/PDESolvers.cpp`](file:///d:/CODING/CM%20PBL/src/unit4/PDESolvers.cpp) |

---

## 3. Architecture & Directory Layout

```text
EES/
├── CMakeLists.txt                # Unified multi-target CMake 3.20+ build system
├── README.md                     # Master project overview & syllabus mapping
├── docs/
│   ├── ARCHITECTURE.md           # Software architecture and design decisions
│   ├── ALGORITHMS.md             # Equations, pseudocode, convergence proofs & complexity
│   ├── MATHEMATICAL_FOUNDATION.md# Analysis theory (IEEE-754, Scarborough, CFL stability)
│   ├── TESTING.md                # Test benchmarks, ground truths & edge cases
│   ├── USER_GUIDE.md             # Complete operational manual & tutorial
│   └── SYLLABUS_MAPPING.md       # Exact 100% topic-by-topic curriculum mapping
├── include/ees/
│   ├── core/                     # Vector, Matrix, ExpressionParser, ErrorAnalysis, Types
│   ├── unit1/                    # Foundations, Root Finding, Optimization
│   ├── unit2/                    # Difference Tables, Interpolation, Integration Quadrature
│   ├── unit3/                    # Linear Systems, LU/Cholesky, Power Method, Splines
│   ├── unit4/                    # ODE Solvers (Euler, RK4, Multi-step) & PDE Solvers
│   ├── comparison/               # Method Benchmark & Comparison Engine
│   ├── advisor/                  # Numerical Method Advisor Heuristic Engine
│   ├── problems/                 # Curated Engineering Problem Catalog
│   └── ui/                       # Visualizer (ASCII plots, heatmaps, HTML exporter) & Workbench
├── src/                          # Implementations matching include/ees/
├── tests/
│   ├── test_core.cpp             # Core math tests
│   ├── test_unit1.cpp            # Unit I automated test suite
│   ├── test_unit2.cpp            # Unit II automated test suite
│   ├── test_unit3.cpp            # Unit III automated test suite
│   └── test_unit4.cpp            # Unit IV automated test suite
└── build/                        # Compiled binaries & test runners
```

---

## 4. Build & Run Instructions

### Prerequisites
- Windows 10 or 11
- CMake 3.20+
- Modern C++ compiler: GCC / MinGW-w64 (GCC 10+) or MSVC (Visual Studio 2019/2022)
- Ninja or Make build tool

### Build Commands (via Ninja)
```powershell
# 1. Configure build directory
cmake -B build -G Ninja

# 2. Compile all targets
ninja -C build

# 3. Run all automated unit tests (100% pass rate)
ctest --test-dir build --output-on-failure
```

### Launching the Application
```powershell
# 1. Direct Command-Line Interface (CLI mode)
.\build\ees.exe root "x^3 - 4*x - 9" --method newton --x0 2.5 --plot --step
.\build\ees.exe root "x^3 - 4*x - 9" --compare
.\build\ees.exe optimize "(x - 3)^2 + 4" --a 0 --b 5 --method golden
.\build\ees.exe integrate "4 / (1 + x^2)" --from 0 --to 1 --method gauss5
.\build\ees.exe solve-system --matrix "[[2,1,-1],[-3,-1,2],[-2,1,2]]" --rhs "[8,-11,-3]" --steps
.\build\ees.exe factorize --matrix "[[4,12],[12,45]]" --method cholesky
.\build\ees.exe ode "y" --x0 0 --y0 1 --xEnd 1 --h 0.1 --method rk4
.\build\ees.exe pde heat --alpha 0.05 --L 1.0 --tMax 0.5

# 2. Interactive Graphical / CLI Workbench
.\build\ees.exe interactive
# or run without arguments:
.\build\ees.exe

# 3. One-Click PBL Viva Demonstration Mode
.\build\ees.exe demo
```

---

## 5. Automated Test Results

All 5 test suites pass with 100% success rate:
```text
Test project D:/CODING/CM PBL/build
    Start 1: CoreMathTests
1/5 Test #1: CoreMathTests ....................   Passed    0.01 sec
    Start 2: Unit1Tests
2/5 Test #2: Unit1Tests .......................   Passed    0.01 sec
    Start 3: Unit2Tests
3/5 Test #3: Unit2Tests .......................   Passed    0.01 sec
    Start 4: Unit3Tests
4/5 Test #4: Unit3Tests .......................   Passed    0.01 sec
    Start 5: Unit4Tests
5/5 Test #5: Unit4Tests .......................   Passed    0.01 sec

100% tests passed, 0 tests failed out of 5
Total Test time (real) = 0.08 sec
```

---

## 6. One-Click PBL Viva Demonstration Showcase

Running `.\build\ees_app.exe --demo` executes the 9 milestone benchmarks sequentially:

```text
========================================================================================
   EES — ENGINEERING EQUATION SOLVER
   ALL-IN-ONE PBL VIVA DEMONSTRATION SHOWCASE
========================================================================================
Executing comprehensive demonstration across all 9 syllabus milestones...

[1/9] ROOT FINDING (Newton-Raphson vs Bisection):
  Root of x^3 - 4x - 9 = 0: 2.706528 in 4 iterations (Order: Quadratic)
[2/9] MULTIVARIATE OPTIMIZATION (Nelder-Mead Simplex):
  Minimum found at: [2.0000, 3.0000] in 39 cycles
[3/9] INTERPOLATION (Gregory-Newton Forward):
  Interpolated y(2.5): 11.625 (Exact: 11.625)
[4/9] NUMERICAL INTEGRATION (Gauss-Legendre 5-Point vs Romberg):
  Computed pi integral: 3.14159264 (Abs Error: 1.3705040214e-08)
[5/9] LINEAR SYSTEM (Gauss Elimination with Partial Pivoting):
  Solution x: [2.0000, 3.0000, -1.0000] (Det: -1.0000000000e+00)
[6/9] EIGENVALUES (Power Method):
  Dominant Eigenvalue: 4.000000 (Vector: [0.7071, 0.7071])
[7/9] NATURAL CUBIC SPLINES (C2 Smoothness):
  Spline interpolation at knot midpoint 1.5: 3.150000
[8/9] ODE SOLUTION (Classic RK4):
  RK4 Final Value: 2.718280 (Global Error: 2.084324e-06)
[9/9] PDE SOLUTION (1D Heat Equation Diffusion):
  Grid Generated: 21 x 11 [CFL r = 2.500000e-01]
PDE Solution Intensity Field [Min=0.000, Max=1.000]:
+----------------------------------------+
|    ::::++++####%%%%@@@@%%%%####++++::: |
|    ::::++++****%%%%%%%%%%%%****++++::: |
|    ::::====****####%%%%####****====::: |
|    ::::====****############****====::: |
|    ::::====++++############++++====::: |
|    ::::====++++************++++====::: |
|    ::::----++++************++++----::: |
|    ....----++++++++****++++++++----... |
|    ....----====++++++++++++====----... |
|    ....----====++++++++++++====----... |
+----------------------------------------+

All 9 syllabus milestones verified and demonstrated successfully.
```

---

## 7. License & Academic Integrity

This software is developed strictly for educational and scientific research as part of the B.Tech Computational Methods Project-Based Learning (PBL) curriculum. Every numerical method and data structure is original C++ code developed from first principles without third-party mathematical solver libraries.

