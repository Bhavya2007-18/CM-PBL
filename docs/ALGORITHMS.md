# ALGORITHMS & MATHEMATICAL IMPLEMENTATION DETAILS

This document presents the mathematical derivations, pseudocode, asymptotic convergence orders, stability proofs, and computational complexities for every algorithm implemented in the **EES (Engineering Equation Solver)** engine.

---

## UNIT I: NUMERICAL FOUNDATIONS & OPTIMIZATION

### 1. Taylor Series Expansion
- **Mathematical Formula**:
  \[
  f(x) = \sum_{k=0}^n \frac{f^{(k)}(x_0)}{k!} (x - x_0)^k + R_n(x)
  \]
  where the Lagrange Remainder is:
  \[
  R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} (x - x_0)^{n+1}, \quad \xi \in (x_0, x)
  \]
- **Implementation**: `ees::unit1::FoundationsLab::computeTaylorSeries`
- **Complexity**: $O(n^2)$ using finite difference higher-order derivative stencils.

### 2. Catastrophic Cancellation & Loss of Significance
- **Problem**: When evaluating $f(x) = \sqrt{x+1} - \sqrt{x}$ for large $x$, significant bits in the mantissa subtract to zero, magnifying round-off error.
- **Stable Form**:
  \[
  \sqrt{x+1} - \sqrt{x} = \frac{1}{\sqrt{x+1} + \sqrt{x}}
  \]
- **Implementation**: `ees::core::ErrorAnalysis::cancellationDirect` vs `cancellationStable`.

### 3. Bisection Method
- **Recurrence**:
  \[
  c_k = \frac{a_k + b_k}{2}
  \]
  If $f(a_k) \cdot f(c_k) < 0 \implies b_{k+1} = c_k$, else $a_{k+1} = c_k$.
- **Error Bound**: $|r - c_k| \le \frac{b - a}{2^k}$.
- **Convergence**: Strictly **Linear** (Order 1, asymptotic error constant $C = 0.5$).
- **Complexity**: $O(\log_2((b-a)/\epsilon))$ iterations.

### 4. Newton-Raphson Method
- **Recurrence**:
  \[
  x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}
  \]
- **Convergence**: **Quadratic** (Order 2):
  \[
  e_{k+1} \approx \frac{f''(r)}{2 f'(r)} e_k^2
  \]
  Doubles the number of accurate decimal places each iteration near a simple root.
- **Failure Modes**: Local extrema ($f'(x) \approx 0$), cycle oscillation, divergence on inflection points.

### 5. Secant Method
- **Recurrence**:
  \[
  x_{k+1} = x_k - f(x_k) \frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}
  \]
- **Convergence**: **Superlinear** (Order $p = \frac{1 + \sqrt{5}}{2} \approx 1.618$).
- **Advantage**: Derivative-free while retaining superlinear speed.

### 6. Nelder-Mead Simplex Algorithm
- **Direct Search Operations**:
  1. **Reflection**: $x_r = \bar{x} + \alpha (\bar{x} - x_{n+1})$ ($\alpha = 1$)
  2. **Expansion**: $x_e = \bar{x} + \gamma (x_r - \bar{x})$ ($\gamma = 2$)
  3. **Contraction**: $x_c = \bar{x} + \rho (x_{n+1} - \bar{x})$ ($\rho = 0.5$)
  4. **Shrink**: $x_i = x_1 + \sigma (x_i - x_1)$ ($\sigma = 0.5$)
- **Complexity**: Gradient-free heuristic, $O(n^2)$ per cycle.

---

## UNIT II: INTERPOLATION & QUADRATURE

### 1. Gregory-Newton Forward & Backward
- **Forward Formula**:
  \[
  P(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \dots, \quad u = \frac{x - x_0}{h}
  \]
- **Backward Formula**:
  \[
  P(x) = y_n + v \nabla y_n + \frac{v(v+1)}{2!} \nabla^2 y_n + \dots, \quad v = \frac{x - x_n}{h}
  \]

### 2. Lagrange & Newton Divided Difference
- **Lagrange Cardinal Polynomial**:
  \[
  P_n(x) = \sum_{i=0}^n y_i \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}
  \]
- **Newton Divided Difference**:
  \[
  P_n(x) = f[x_0] + \sum_{k=1}^n f[x_0, \dots, x_k] \prod_{j=0}^{k-1} (x - x_j)
  \]

### 3. Quadrature Formulas & Errors
- **Trapezoidal**: $\frac{h}{2} [f_0 + 2\sum f_i + f_n]$, Truncation Error $E_T = -\frac{(b-a)h^2}{12} f''(\xi) = O(h^2)$
- **Simpson's 1/3**: $\frac{h}{3} [f_0 + 4\sum_{odd} f_i + 2\sum_{even} f_i + f_n]$, Truncation Error $E_S = -\frac{(b-a)h^4}{180} f^{(4)}(\xi) = O(h^4)$
- **Simpson's 3/8**: $\frac{3h}{8} [f_0 + 3f_1 + 3f_2 + 2f_3 + \dots + f_n]$, Error $O(h^4)$
- **Romberg Extrapolation**:
  \[
  R(k, j) = R(k, j-1) + \frac{R(k, j-1) - R(k-1, j-1)}{4^j - 1}
  \]
- **Gauss-Legendre Quadrature**:
  \[
  \int_a^b f(x) dx \approx \frac{b-a}{2} \sum_{i=1}^n w_i f\left( \frac{b-a}{2} t_i + \frac{a+b}{2} \right)
  \]
  Integrates polynomials up to degree $2n - 1$ exactly using only $n$ points.

---

## UNIT III: LINEAR SYSTEMS, FACTORIZATION & SPLINES

### 1. Gauss Elimination with Partial Pivoting
- **Row Pivot**: $\max_{k \le i \le n} |A_{ik}|$ swapped to $(k, k)$.
- **Multiplier**: $m_{ik} = A_{ik} / A_{kk}$, row update $R_i \leftarrow R_i - m_{ik} R_k$.
- **Complexity**: $\frac{2}{3} n^3$ additions/multiplications.

### 2. Matrix Factorizations
- **Doolittle**: $A = LU$ where $L_{ii} = 1$.
- **Crout**: $A = LU$ where $U_{ii} = 1$.
- **Cholesky**: $A = L L^T$ where $L_{ii} = \sqrt{A_{ii} - \sum_{k=1}^{i-1} L_{ik}^2}$. Complexity $\frac{1}{3} n^3$ FLOPs.

### 3. Natural Cubic Splines ($C^2$ Continuous)
- Piecewise cubic: $S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3$.
- Natural boundary: $S''(x_0) = 0, S''(x_n) = 0$.
- Solves tridiagonal system in $O(n)$ time using the Thomas Algorithm.

---

## UNIT IV: ODES & PDES

### 1. Runge-Kutta 4th-Order (RK4)
- **Formulas**:
  \[
  k_1 = f(x_n, y_n)
  \]
  \[
  k_2 = f\left(x_n + \frac{h}{2}, y_n + \frac{h}{2} k_1\right)
  \]
  \[
  k_3 = f\left(x_n + \frac{h}{2}, y_n + \frac{h}{2} k_2\right)
  \]
  \[
  k_4 = f(x_n + h, y_n + h k_3)
  \]
  \[
  y_{n+1} = y_n + \frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4)
  \]
- **Global Error**: $O(h^4)$, Local Error $O(h^5)$.

### 2. Predictor-Corrector Methods
- **Adams-Bashforth-Moulton**:
  - Predictor (AB4): $y_{n+1}^* = y_n + \frac{h}{24} [55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3}]$
  - Corrector (AM3): $y_{n+1} = y_n + \frac{h}{24} [9 f_{n+1} + 19 f_n - 5 f_{n-1} + f_{n-2}]$
- **Milne-Simpson**:
  - Predictor: $y_{n+1}^* = y_{n-3} + \frac{4h}{3} [2 f_n - f_{n-1} + 2 f_{n-2}]$
  - Corrector: $y_{n+1} = y_{n-1} + \frac{h}{3} [f_{n-1} + 4 f_n + f_{n+1}]$

### 3. PDE Solvers
- **Heat Equation (Parabolic)**:
  \[
  u_i^{n+1} = u_i^n + r (u_{i+1}^n - 2u_i^n + u_{i-1}^n), \quad r = \frac{\alpha \Delta t}{\Delta x^2} \le \frac{1}{2}
  \]
- **Laplace Equation (Elliptic)**:
  \[
  u_{i,j}^{(k+1)} = (1 - \omega) u_{i,j}^{(k)} + \frac{\omega}{4} [u_{i+1,j}^{(k)} + u_{i-1,j}^{(k+1)} + u_{i,j+1}^{(k)} + u_{i,j-1}^{(k+1)}]
  \]
- **Wave Equation (Hyperbolic)**:
  \[
  u_i^{n+1} = 2(1 - C^2) u_i^n + C^2 (u_{i+1}^n + u_{i-1}^n) - u_i^{n-1}, \quad C = \frac{c \Delta t}{\Delta x} \le 1
  \]

