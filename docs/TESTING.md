# TESTING & VALIDATION STRATEGY

---

## 1. Test Suite Architecture

The test infrastructure of **EES (Engineering Equation Solver)** is organized into modular automated test suites integrated directly with CTest and CMake:

| Test Target | Source File | Scope of Verification | Pass Status |
|---|---|---|---|
| `CoreMathTests` | `tests/test_core.cpp` | `Vector`, `Matrix`, `ExpressionParser`, `ErrorAnalysis`, `NumericalResult` | **PASSED (5/5)** |
| `Unit1Tests` | `tests/test_unit1.cpp` | Taylor Series, Rolle's, MVT, Bisection, Newton, Secant, Fibonacci, Golden Section, Nelder-Mead | **PASSED (3/3)** |
| `Unit2Tests` | `tests/test_unit2.cpp` | Difference Tables, Gregory-Newton F/B, Lagrange, Divided Difference, Quadratures (Trapezoidal, Simpson 1/3, 3/8, Romberg, Gauss 2-5 pt) | **PASSED (2/2)** |
| `Unit3Tests` | `tests/test_unit3.cpp` | Gauss Elimination with Pivoting, Gauss-Jordan RREF, Matrix Inverse, Doolittle LU, Crout LU, Cholesky $LL^T$, Power Method, Splines | **PASSED (4/4)** |
| `Unit4Tests` | `tests/test_unit4.cpp` | ODEs (Euler, Heun, RK4, Adams-Bashforth, Milne), PDEs (Heat, Laplace SOR, Wave Equation) | **PASSED (2/2)** |

Total test pass rate: **100% (16/16 test suites and sub-verifications)**.

---

## 2. Benchmark Problems & Ground Truth

### Root Finding Benchmarks
- Equation: $f(x) = x^3 - 4x - 9 = 0$.
  - Interval: $[2, 3]$.
  - Ground Truth Root: $r \approx 2.70652796$.
  - Tolerance: $\epsilon = 10^{-7}$.
  - Verification: Bisection converges linearly in 18 iterations; Newton-Raphson converges quadratically in 4 iterations; Secant converges superlinearly in 6 iterations.

### Quadrature Benchmarks
- Integral: $\int_0^1 \frac{4}{1 + x^2} dx = \pi \approx 3.141592653589793$.
  - 100-step Trapezoidal: Absolute error $< 10^{-4}$.
  - 100-step Simpson's 1/3: Absolute error $< 10^{-7}$.
  - 6th-order Romberg: Absolute error $< 10^{-9}$.
  - 5-point Gauss-Legendre: Absolute error $< 10^{-7}$ using only 5 evaluations.

### Linear Systems & Matrix Factorizations
- System:
  \[
  \begin{bmatrix} 2 & 1 & -1 \\ -3 & -1 & 2 \\ -2 & 1 & 2 \end{bmatrix} \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{bmatrix} 8 \\ -11 \\ -3 \end{bmatrix}
  \]
  - Ground Truth Solution: $x = 2.0, y = 3.0, z = -1.0$.
  - Determinant: $-1.0$.
  - Matrix Inversion: Verified $A \times A^{-1} = I$ to within $\|A A^{-1} - I\|_\infty < 10^{-10}$.
- Cholesky Factorization on Symmetric Positive Definite Matrix:
  \[
  A = \begin{bmatrix} 4 & 12 \\ 12 & 45 \end{bmatrix} \implies L = \begin{bmatrix} 2 & 0 \\ 6 & 3 \end{bmatrix}, \quad \|A - L L^T\|_F < 10^{-12}.
  \]

### Ordinary Differential Equations
- Initial Value Problem: $\frac{dy}{dx} = y, \quad y(0) = 1, \quad x \in [0, 1]$.
  - Exact Analytic Solution: $y(x) = e^x \implies y(1) \approx 2.718281828$.
  - Step size $h = 0.1$:
    - Euler: Final error $\approx 0.12$.
    - Heun (Euler Predictor-Corrector): Final error $\approx 0.006$.
    - Runge-Kutta 4th-Order (RK4): Final error $\approx 2.08 \times 10^{-6}$.
    - Adams-Bashforth-Moulton 4-step: Final error $\approx 1.5 \times 10^{-5}$.
    - Milne's Predictor-Corrector: Final error $\approx 1.2 \times 10^{-5}$.

### Partial Differential Equations
- 1D Heat Conduction: $\alpha = 0.1, L = 1.0, u(x, 0) = \sin(\pi x)$, Dirichlet boundaries $0$.
  - Verification: Center temperature starts at $1.0$ and diffuses monotonically without numerical oscillations when $r = \frac{\alpha \Delta t}{\Delta x^2} = 0.25 \le 0.5$.
- 2D Laplace Equation: $15 \times 15$ grid, Top $= 100^\circ\text{C}$, others $0^\circ\text{C}$.
  - Verification: Solved via SOR ($\omega = 1.3$) in under 60 cycles; satisfies discrete maximum principle $0 \le u_{i,j} \le 100$.

---

## 3. Boundary & Edge Case Handling

1. **Bracketing Violation**: Passing an interval where $f(a) \cdot f(b) > 0$ to Bisection returns `ConvergenceStatus::InvalidInput` with a clear explanation rather than an infinite loop.
2. **Vanishing Derivative**: When $f'(x) \to 0$ in Newton-Raphson, division by zero is caught safely and reports `ConvergenceStatus::DivisionByZero`.
3. **Singular Matrix**: Zero pivots in Gaussian elimination trigger partial pivoting row swaps; if all column candidates are zero, reports `ConvergenceStatus::SingularMatrix`.
4. **CFL Instability**: In the Heat Equation, if the user supplies $r > 0.5$, EES detects $r > 0.5$ and flags `[UNSTABLE - Exceeds CFL limit!]`.

---

## 4. Running the Tests

To run the automated suite from PowerShell or Command Prompt:

```powershell
# Run individual test executables
.\build\test_core.exe
.\build\test_unit1.exe
.\build\test_unit2.exe
.\build\test_unit3.exe
.\build\test_unit4.exe

# Or run the entire CTest harness
ctest --test-dir build --output-on-failure
```

