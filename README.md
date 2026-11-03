# EES — Engineering Equation Solver
### Numerical Methods & Engineering Analysis Workbench
**A Unified Computational Methods / Numerical Methods PBL Platform (C++17)**

![C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Build-Tool](https://img.shields.io/badge/Build-CMake%20%7C%20Ninja-success.svg)
![Tests](https://img.shields.io/badge/Tests-16%2F16%20Passing%20(100%25)-brightgreen.svg)
![Curriculum](https://img.shields.io/badge/Curriculum-Units%20I--IV%20(100%25)-orange.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%2010%2F11-blueviolet.svg)

---

## 1. Project Overview & Philosophy

> **Solve → Explain → Visualize → Compare → Understand**

**EES (Engineering Equation Solver)** is an engineering-grade scientific computing application engineered from first principles in modern C++ (C++17/20). Rather than a collection of disconnected practical scripts, EES provides an integrated **Numerical Methods Engineering Workbench** designed for university Computational Methods Project-Based Learning (PBL) demonstrations, viva defense, and engineering design.

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

## 3. Directory Layout

```text
EES/
├── CMakeLists.txt                # Unified multi-target CMake 3.20+ build configuration
├── README.md                     # Master project documentation (this file)
├── LICENSE                       # MIT Open Source License
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
│   └── ui/                       # Visualizer (ASCII plots, heatmaps, HTML exporter), Workbench & CLI
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

## 4. How to Build & Run

### Prerequisites
- **Operating System**: Windows 10 or Windows 11 (64-bit)
- **Compiler**: GCC / MinGW-w64 (GCC 10+) or MSVC (Visual Studio 2019/2022) or Clang
- **Build System**: CMake 3.20+ and Ninja (or Make)

### Step 1: Configure & Compile

Open **PowerShell** or **Command Prompt** in the project directory (`d:\CODING\CM PBL`):

```powershell
# 1. Generate Ninja build files with CMake
cmake -B build -G Ninja

# 2. Compile the entire application and test suite
ninja -C build
```

*(Alternatively, if using Visual Studio generator: `cmake -B build` and `cmake --build build --config Release`)*

---

### Step 2: How to Run the Application

The build produces two executable targets:
- `build/ees.exe` — The primary Unified CLI and Workbench binary
- `build/ees_app.exe` — Alternative binary name

There are **three distinct ways** to run EES:

#### Mode A: Interactive Desktop Workbench (Recommended for Exploring)
Run the executable with no arguments (or with `interactive`):

```powershell
.\build\ees.exe
# OR
.\build\ees.exe interactive
```

This launches the terminal scientific workbench with full menu navigation:
```text
========================================================================================
   EES — ENGINEERING EQUATION SOLVER
   MAIN WORKBENCH DASHBOARD
========================================================================================
  [1] UNIT I:  Foundations, Error Analysis & Computer Arithmetic Lab
  [2] UNIT I:  Root-Finding Laboratory (Bisection, Newton, Secant)
  [3] UNIT I:  Optimization (Fibonacci, Golden Section, Nelder-Mead)
  [4] UNIT II: Interpolation & Finite Differences Laboratory
  [5] UNIT II: Numerical Integration & Quadrature Laboratory
  [6] UNIT III: Linear Systems, Pivoting, LU & Cholesky Decomposition
  [7] UNIT III: Dominant Eigenvalues (Power Method)
  [8] UNIT III: Spline Approximation (Linear, Quadratic, Natural Cubic, B-Spline)
  [9] UNIT IV: Numerical Solution of ODEs (Euler, RK4, Predictor-Corrector)
 [10] UNIT IV: Numerical Solution of PDEs (Heat, Wave, Laplace/Poisson)
 [11] METHOD COMPARISON BENCHMARK ENGINE
 [12] ENGINEERING PROBLEM CATALOG (Mechanical, Civil, Electrical, Thermal)
 [13] NUMERICAL METHOD ADVISOR (Intelligent Recommendation Engine)
 [14] ALL-IN-ONE PBL VIVA DEMONSTRATION MODE (One-Click Showcase)
  [0] Exit EES
```

#### Mode B: One-Click PBL Viva Demonstration (Best for Viva Defense)
To immediately execute all 9 core syllabus milestone benchmarks sequentially without pausing:

```powershell
.\build\ees.exe demo
# OR
.\build\ees.exe --demo
```

#### Mode C: Direct Command-Line Interface (CLI Mode)
You can directly solve equations, integrate, invert matrices, solve differential equations, or compare methods from your terminal:

```powershell
# 1. Root Finding: solve custom equation with quadratic Newton-Raphson & plot curve
.\build\ees.exe root "x^3 - 4*x - 9" --method newton --x0 2.5 --plot --step

# 2. Compare Bisection vs Newton vs Secant on the same equation
.\build\ees.exe root "x^3 - 4*x - 9" --compare

# 3. 1D Optimization (Golden Section Search)
.\build\ees.exe optimize "(x - 3)^2 + 4" --a 0 --b 5 --method golden

# 4. Numerical Integration (5-point Gauss-Legendre Quadrature)
.\build\ees.exe integrate "4 / (1 + x^2)" --from 0 --to 1 --method gauss5 --exact 3.141592653589793

# 5. Linear Systems: Solve Ax = b and inspect step-by-step row operations
.\build\ees.exe solve-system --matrix "[[2,1,-1],[-3,-1,2],[-2,1,2]]" --rhs "[8,-11,-3]" --steps

# 6. Matrix Factorization (Cholesky A = L*L^T on symmetric positive-definite matrix)
.\build\ees.exe factorize --matrix "[[4,12],[12,45]]" --method cholesky

# 7. Dominant Eigenvalue via Power Method with Rayleigh Quotient
.\build\ees.exe eigen --matrix "[[2,2],[1,3]]" --plot

# 8. Natural Cubic Spline Interpolation
.\build\ees.exe spline --x "0,1,2,3" --y "0,1,8,27" --eval 1.5

# 9. Initial Value ODE: Solve dy/dx = y using Runge-Kutta 4th-Order (RK4)
.\build\ees.exe ode "y" --x0 0 --y0 1 --xEnd 1 --h 0.1 --method rk4 --table

# 10. Compare all ODE solvers (Euler, Heun, RK4, Adams-Bashforth, Milne)
.\build\ees.exe ode "y" --x0 0 --y0 1 --xEnd 1 --h 0.1 --compare

# 11. 1D Parabolic Heat Equation PDE (with CFL stability check & ASCII heatmap)
.\build\ees.exe pde heat --alpha 0.05 --L 1.0 --tMax 0.5

# 12. Floating-Point Stability & Catastrophic Cancellation Demonstration
.\build\ees.exe arithmetic

# 13. Heuristic Numerical Method Advisor
.\build\ees.exe advise root --derivative --bracket

# 14. Display Full CLI Help & Flags
.\build\ees.exe --help
```

---

### Step 3: Running Automated Test Suites

EES includes 5 comprehensive test executables integrated with **CTest**:

```powershell
# Run the complete test suite via CTest
ctest --test-dir build --output-on-failure

# Or run individual test suites directly:
.\build\test_core.exe     # Vector, Matrix, Parser, ErrorAnalysis
.\build\test_unit1.exe    # Taylor, Roots, Optimization
.\build\test_unit2.exe    # Differences, Interpolation, Quadratures
.\build\test_unit3.exe    # Linear Systems, LU, Cholesky, Splines
.\build\test_unit4.exe    # ODEs (Euler, RK4, Milne) & PDEs (Heat, Wave, Laplace)
```

**Test Verification Summary**:
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
Total Test time (real) = 0.06 sec
```

---

## 5. Sample Visualizations in Terminal

### Terminal 2D Function Curve with Root Marker `[O]`:
```text
=== Root Finding Plot ===
  y_max = 21.4726
  |        |                                                       *|
  |        |                                                      * |
  |        |                                                    **  |
  |        |                                                   *    |
  |        |                                                 **     |
  |        |                                                *       |
  |        |                                              **        |
  |        |                                            **          |
  |--------+------------------------------------------O*------------|
  |        |                                       ***              |
  |        |                                    ***                 |
  |******* |                                ****                    |
  |       ***********                 ******                        |
  |        |         *****************                              |
  |        |                                                        |
  y_min = -13.6737
  x: [-0.5000 ... 3.5000]
  Legend: [*] f(x) curve   [O] Root / Iteration Marker   [+] Origin
```

### Semi-Log Convergence Trajectory Plot:
```text
=== Convergence Trajectory (Iteration vs Error) ===
  log10(Error) [Top = 10^-0.6]
  |#                                                                |
  |                     #                                           |
  |                                                                 |
  |                                          #                      |
  |                                                                 |
  |                                                                #|
  log10(Error) [Bottom = 10^-8]
  Iteration: 1 -> 4
```

### 1D Heat Equation Spatial-Temporal Heatmap:
```text
=== PDE Solution (Parabolic: 1D Heat Equation) ===
Numerical Method: FTCS Explicit Finite Difference Scheme
Grid Dimensions:  51 x 21
CFL Stability Number: 0.2000 [STABLE]

PDE Solution Intensity Field [Min=0.000, Max=100.000]:
+---------------------------------------------+
|   ..::==++***##%%%%%%@@@%%%%%%##***++==::.. |
|   ..::--++*****##%%%%%%%%%%%##*****++--::.. |
|   ..::--==+++**####%%%%%%%####**+++==--::.. |
|   ..::--==+++**###############**+++==--::.. |
|   ..::--==+++****###########****+++==--::.. |
|   ..::--==+++++****#######****+++++==--::.. |
|   ..::--=====++***************++=====--::.. |
|   ..::----===++++***********++++===----::.. |
|     ..::--===++++++*******++++++===--::..   |
|     ..::--=====++++++***++++++=====--::..   |
|     ..::--=====+++++++++++++++=====--::..   |
|     ..::-----====+++++++++++====-----::..   |
+---------------------------------------------+
```

---

## 6. Curated Engineering Problem Catalog

EES provides pre-built engineering benchmarks ready to solve via menu Option 12:

1. **Mechanical Engineering (`MECH-01`)**:
   - **Problem**: Van der Waals equation of state for real gases: $(P + a/V^2)(V - b) = RT$.
   - **Application**: Carbon dioxide molar volume at $50\text{ atm}$, $320\text{ K}$.
   - **Method**: Newton-Raphson quadratic solve ($V \approx 0.000378\text{ m}^3/\text{mol}$).
2. **Civil Engineering (`CIVIL-01`)**:
   - **Problem**: 3-bar pin-jointed planar structural truss equilibrium.
   - **Application**: Joint load resolution under $10\text{ kN}$ vertical load.
   - **Method**: Gauss Elimination with Partial Pivoting.
3. **Electrical Engineering (`ELEC-01`)**:
   - **Problem**: Nonlinear Shockley diode equation in series circuit: $(V_s - V_D)/R = I_s(e^{V_D/V_t} - 1)$.
   - **Application**: Forward operating voltage $V_D \approx 0.684\text{ V}$.
   - **Method**: Newton-Raphson with central numerical derivatives.
4. **Thermal Engineering (`THERM-01`)**:
   - **Problem**: Transient heat conduction in aluminum rod quenched at both ends: $\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}$.
   - **Method**: Parabolic finite-difference FTCS scheme ($r \le 0.5$).
5. **General Engineering (`GEN-01`)**:
   - **Problem**: Lubricating oil viscosity vs temperature calibration.
   - **Method**: Natural Cubic Spline ($C^2$ continuous, minimal bending energy).

---

## 7. Documentation Index

For exhaustive theoretical derivations, pseudocode, and usage guides:
- [`docs/SYLLABUS_MAPPING.md`](file:///d:/CODING/CM%20PBL/docs/SYLLABUS_MAPPING.md) — 100% Curriculum syllabus mapping table.
- [`docs/ARCHITECTURE.md`](file:///d:/CODING/CM%20PBL/docs/ARCHITECTURE.md) — Software architecture, layering, and design patterns.
- [`docs/ALGORITHMS.md`](file:///d:/CODING/CM%20PBL/docs/ALGORITHMS.md) — Mathematical formulas, convergence orders, and complexities ($O(n^3)$, $O(h^4)$).
- [`docs/MATHEMATICAL_FOUNDATION.md`](file:///d:/CODING/CM%20PBL/docs/MATHEMATICAL_FOUNDATION.md) — Numerical mathematics theory (IEEE-754, Scarborough, CFL stability).
- [`docs/TESTING.md`](file:///d:/CODING/CM%20PBL/docs/TESTING.md) — Test plan, benchmark values, and edge-case handling.
- [`docs/USER_GUIDE.md`](file:///d:/CODING/CM%20PBL/docs/USER_GUIDE.md) — Comprehensive user guide with command reference.

---

## 8. License & Academic Integrity

This software is developed strictly for educational and scientific research as part of the B.Tech Computational Methods Project-Based Learning (PBL) curriculum. Every numerical method and data structure is original C++ code developed from first principles without third-party mathematical solver libraries.

Licensed under the [MIT License](LICENSE).
