# MATHEMATICAL FOUNDATIONS OF COMPUTATIONAL METHODS

---

## 1. Floating-Point Arithmetic and Precision Limits

Modern scientific computation relies upon IEEE-754 binary floating-point representation. A real number is approximated as:
\[
x = (-1)^s \times (1.m)_2 \times 2^{e - \text{bias}}
\]
- Single Precision (`float`, 32-bit): 1 sign bit, 8 exponent bits (bias 127), 23 fraction bits ($\approx 7$ decimal significant digits). Machine epsilon $\epsilon_m \approx 1.19 \times 10^{-7}$.
- Double Precision (`double`, 64-bit): 1 sign bit, 11 exponent bits (bias 1023), 52 fraction bits ($\approx 16$ decimal significant digits). Machine epsilon $\epsilon_m \approx 2.22 \times 10^{-16}$.

### Scarborough Criterion
To guarantee that an approximation $x^*$ is accurate to at least $m$ significant digits:
\[
\text{Relative Error } \epsilon_r = \frac{|x - x^*|}{|x|} \le \frac{1}{2} \times 10^{-m}
\]
Hence, the number of guaranteed significant digits is computed in EES as:
\[
m = -\left\lfloor \log_{10}(2 \epsilon_r) \right\rfloor
\]

---

## 2. Order of Convergence in Root Finding

Let a sequence $\{x_k\}$ converge to root $r$ with error $e_k = |x_k - r|$. If positive constants $C$ and $p$ exist such that:
\[
\lim_{k \to \infty} \frac{e_{k+1}}{e_k^p} = C
\]
then $p$ is the **order of convergence**, and $C$ is the **asymptotic error constant**:
- $p = 1$: **Linear convergence** (Bisection: $C = 0.5$).
- $1 < p < 2$: **Superlinear convergence** (Secant: $p = \frac{1+\sqrt{5}}{2} \approx 1.618$).
- $p = 2$: **Quadratic convergence** (Newton-Raphson: $C = \frac{|f''(r)|}{2|f'(r)|}$).
- $p = 3$: **Cubic convergence** (Halley's method).

---

## 3. Finite Differences & Interpolation Theory

For equispaced data nodes $x_i = x_0 + i h$:
- Forward Difference: $\Delta y_i = y_{i+1} - y_i$.
- Backward Difference: $\nabla y_i = y_i - y_{i-1}$.
- Shift operator $E$: $E y_i = y_{i+1} = (1 + \Delta) y_i \implies E^u = (1 + \Delta)^u$.

By binomial series expansion of $(1 + \Delta)^u$:
\[
y(x_0 + u h) = E^u y_0 = \sum_{k=0}^\infty \binom{u}{k} \Delta^k y_0 = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \dots
\]
This derives the **Gregory-Newton Forward** formula from operator algebra.

### Polynomial Error Remainder
If $f \in C^{n+1}[a, b]$, the interpolation error at any point $x$ is:
\[
f(x) - P_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{i=0}^n (x - x_i), \quad \xi \in (a, b)
\]
This explicitly shows why high-degree polynomial interpolation on equispaced grids diverges at boundaries (**Runge's Phenomenon**), justifying the transition to piecewise **Natural Cubic Splines** and **B-Splines**.

---

## 4. Stability in Numerical PDEs

The 1D Heat Equation $\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}$ discretized using the Forward-Time Central-Space (FTCS) finite-difference scheme yields the amplified matrix equation:
\[
u^{n+1} = G u^n
\]
By Von Neumann stability analysis, replacing $u_j^n$ by $e^{i k j \Delta x}$:
\[
\xi = 1 - 4 r \sin^2\left(\frac{k \Delta x}{2}\right)
\]
Stability requires $|\xi| \le 1$ for all wave numbers $k$:
\[
-1 \le 1 - 4 r \sin^2\left(\frac{k \Delta x}{2}\right) \le 1 \implies 4 r \le 2 \implies r = \frac{\alpha \Delta t}{\Delta x^2} \le \frac{1}{2}
\]
Violating this Courant-Friedrichs-Lewy (CFL) limit causes instantaneous numerical explosion (exponential growth of round-off errors).

