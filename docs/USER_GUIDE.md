# EES USER MANUAL & OPERATIONAL GUIDE

---

## 1. Quickstart & Command Line Interface

The **EES** application provides both a powerful direct **Command-Line Interface (CLI)** and an interactive **Desktop Workbench Shell**.

### 1.1 Direct CLI Mode

Run individual algorithms directly with custom equations and flags:

```powershell
# Root Finding
.\build\ees.exe root "x^3 - 4*x - 9" --method newton --x0 2.5 --plot --step
.\build\ees.exe root "x^3 - 4*x - 9" --compare

# 1D Optimization
.\build\ees.exe optimize "(x - 3)^2 + 4" --a 0 --b 5 --method golden

# Numerical Integration
.\build\ees.exe integrate "4 / (1 + x^2)" --from 0 --to 1 --method gauss5 --exact 3.141592653589793

# Interpolation
.\build\ees.exe interpolate --x "1,2,3,4" --y "0,5,22,57" --eval 2.5 --method newton-forward --diff-table

# Linear Systems & Factorization
.\build\ees.exe solve-system --matrix "[[2,1,-1],[-3,-1,2],[-2,1,2]]" --rhs "[8,-11,-3]" --steps
.\build\ees.exe factorize --matrix "[[4,12],[12,45]]" --method cholesky
.\build\ees.exe eigen --matrix "[[2,2],[1,3]]" --plot

# Differential Equations
.\build\ees.exe ode "y" --x0 0 --y0 1 --xEnd 1 --h 0.1 --method rk4
.\build\ees.exe pde heat --alpha 0.05 --L 1.0 --tMax 0.5

# Method Advisor
.\build\ees.exe advise root --derivative --bracket

# Floating-Point Stability & Loss of Significance Lab
.\build\ees.exe arithmetic
```

### 1.2 Interactive Workbench & Demo Mode

```powershell
# Launch the interactive scientific workbench
.\build\ees.exe interactive
# or simply run without arguments:
.\build\ees.exe

# Run the 9-milestone PBL Viva Demonstration showcase non-interactively
.\build\ees.exe demo
```

---

## 2. Main Workbench Navigation

Upon startup, the desktop-style scientific dashboard appears:

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

---

## 3. Module Walkthroughs

### Option 1: Numerical Foundations & Computer Arithmetic
- **Computer Arithmetic Lab**: Demonstrates floating-point bit decomposition (sign, exponent, mantissa), IEEE-754 machine epsilon, Scarborough significant digit criterion, and catastrophic cancellation when subtracting nearly equal numbers (e.g. $\sqrt{x+1} - \sqrt{x}$ vs $\frac{1}{\sqrt{x+1}+\sqrt{x}}$).
- **Taylor Series**: Enter any function $f(x)$ (e.g. `exp(x)`, `sin(x)`), expansion point $x_0$, order $n$, and evaluation point $x$. Displays term-by-term convergence and error table.
- **Rolle's & Mean Value Theorems**: Numerically searches for and displays interior points $c \in (a, b)$ where $f'(c) = 0$ (Rolle's) or $f'(c) = \frac{f(b)-f(a)}{b-a}$ (MVT).

### Option 2: Root-Finding Laboratory
- Solves arbitrary user equations $f(x) = 0$ using Bisection, Newton-Raphson, and Secant methods.
- Generates:
  - Iteration table with $x_k, f(x_k)$, absolute error, relative error, and auxiliary parameters ($f'(x)$, interval width).
  - Terminal 2D plot of function curve with root location marked `[O]`.
  - Semi-log convergence plot ($\log_{10}(\text{Error})$ vs Iteration).
  - Search interval shrinkage bar visualization.

### Option 4: Interpolation & Finite Differences
- Computes complete forward, backward, and divided difference tables.
- Evaluates interpolating polynomials using Gregory-Newton Forward/Backward, Lagrange, and Newton Divided Difference.
- Displays full symbolic polynomial formulas.

### Option 5: Numerical Integration
- Solves definite integrals using Trapezoidal, Simpson's 1/3, Simpson's 3/8, Romberg Extrapolation, and Gauss-Legendre Quadrature (2, 3, 4, 5-point).
- Displays full Romberg extrapolation triangular table.

### Option 6: Linear Systems & Matrix Factorization
- Solves systems $Ax = b$ via Gaussian Elimination with Partial Pivoting and Gauss-Jordan (Reduced Row Echelon Form).
- Outputs every row operation step (e.g. `Swap R1 <-> R2`, `R2 <- R2 - (0.5)*R1`).
- Factorizes matrices into Doolittle $LU$, Crout $LU$, and Cholesky $LL^T$ with Frobenius reconstruction error.
- Interactive Pivoting Lab demonstrates catastrophic instability when pivoting is omitted on ill-conditioned systems.

### Option 7: Power Method Dominant Eigenvalues
- Computes dominant eigenvalue $\lambda_1$ and normalized eigenvector $v_1$ using Rayleigh quotient iterations.
- Plots convergence trajectory in real time.

### Option 8: Spline Approximations
- Computes Linear, Quadratic, and Natural Cubic Splines ($C^2$ continuous).
- Generates piecewise cubic polynomial coefficients:
  \[
  S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
  \]
- Evaluates B-Splines using Cox-de Boor recursive blending.

### Option 9 & 10: ODE & PDE Solvers
- **ODEs**: Solves initial value problems using Euler, Heun (Predictor-Corrector), RK4, Adams-Bashforth 4-step, and Milne's method.
- **PDEs**: Solves 1D Heat Equation (parabolic), 2D Laplace Equation (elliptic), and 1D Wave Equation (hyperbolic). Displays high-density ASCII/Unicode heatmaps directly in the terminal!

### Option 11: Method Comparison Benchmark Engine
- Solves the same problem simultaneously with multiple numerical methods.
- Displays a comparative matrix of results, iterations, final error, execution time, and convergence order.

### Option 12: Engineering Problem Catalog
- Pre-formulated real-world engineering benchmarks across Mechanical (Van der Waals CO2 gas volume), Civil (planar truss forces), Electrical (Shockley diode operating point), Thermal (cooling rod transient diffusion), and General Engineering (lubricating oil viscosity).

### Option 13: Numerical Method Advisor
- Interactive expert system that asks for problem characteristics (e.g. interval availability, derivative computability, matrix symmetry) and outputs mathematically justified algorithmic recommendations.

### Option 14: One-Click PBL Viva Demo Mode
- Immediately executes all 9 syllabus milestones in sequence, providing a comprehensive viva demonstration.

