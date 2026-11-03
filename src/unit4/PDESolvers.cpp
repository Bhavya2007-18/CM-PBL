#include "ees/unit4/PDESolvers.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace ees::unit4 {

std::string PDESolution::summary() const {
    std::ostringstream ss;
    ss << "=== PDE Solution (" << pdeType << ") ===\n";
    ss << "Numerical Method: " << method << "\n";
    ss << "Grid Dimensions:  " << grid.rows() << " x " << grid.cols() << "\n";
    if (parameterR > 0.0) {
        ss << "CFL Stability Number: " << std::fixed << std::setprecision(4) << parameterR;
        if (isStable) {
            ss << " [STABLE]\n";
        } else {
            ss << " [UNSTABLE - Exceeds CFL limit!]\n";
        }
    }
    if (iterations > 0) {
        ss << "Iterative Solver Cycles: " << iterations << "\n";
    }
    return ss.str();
}

std::string PDESolution::asciiHeatmap(int displayCols, int displayRows) const {
    if (grid.empty()) return "Empty PDE grid.\n";

    size_t rows = grid.rows();
    size_t cols = grid.cols();

    // Find min and max
    double minVal = grid(0, 0);
    double maxVal = grid(0, 0);
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            minVal = std::min(minVal, grid(r, c));
            maxVal = std::max(maxVal, grid(r, c));
        }
    }

    double range = maxVal - minVal;
    if (std::abs(range) < 1e-12) range = 1.0;

    // Density palette from low to high: " .:-=+*#%@"
    const std::string palette = " .:-=+*#%@";
    const int paletteLen = static_cast<int>(palette.length());

    std::ostringstream ss;
    ss << "PDE Solution Intensity Field [Min=" << std::fixed << std::setprecision(3)
       << minVal << ", Max=" << maxVal << "]:\n";
    ss << "+" << std::string(displayCols, '-') << "+\n";

    for (int dr = 0; dr < displayRows; ++dr) {
        size_t r = (dr * (rows - 1)) / (displayRows > 1 ? (displayRows - 1) : 1);
        ss << "|";
        for (int dc = 0; dc < displayCols; ++dc) {
            size_t c = (dc * (cols - 1)) / (displayCols > 1 ? (displayCols - 1) : 1);
            double val = grid(r, c);
            int idx = static_cast<int>(((val - minVal) / range) * (paletteLen - 1));
            idx = std::clamp(idx, 0, paletteLen - 1);
            ss << palette[idx];
        }
        ss << "|\n";
    }
    ss << "+" << std::string(displayCols, '-') << "+\n";
    return ss.str();
}

PDESolution PDESolverLab::solveHeatEquation1D(
    double alpha,
    double L,
    double tMax,
    size_t numX,
    size_t numT,
    core::RealFunction initialCondition,
    double boundaryLeft,
    double boundaryRight
) {
    PDESolution sol;
    sol.pdeType = "Parabolic: 1D Heat Equation";
    sol.method = "FTCS Explicit Finite Difference Scheme";

    double dx = L / static_cast<double>(numX - 1);
    double dt = tMax / static_cast<double>(numT - 1);
    double r = alpha * dt / (dx * dx);

    sol.parameterR = r;
    sol.isStable = (r <= 0.5);

    sol.grid = core::Matrix::zeros(numT, numX);
    sol.xCoords = core::Vector::linspace(0.0, L, numX);
    sol.yCoordsOrTime = core::Vector::linspace(0.0, tMax, numT);

    // Initial condition at t = 0 (row 0)
    for (size_t i = 0; i < numX; ++i) {
        sol.grid(0, i) = initialCondition(sol.xCoords[i]);
    }
    // Boundary conditions
    sol.grid(0, 0) = boundaryLeft;
    sol.grid(0, numX - 1) = boundaryRight;

    // Time-stepping forward
    for (size_t n = 0; n < numT - 1; ++n) {
        sol.grid(n + 1, 0) = boundaryLeft;
        sol.grid(n + 1, numX - 1) = boundaryRight;

        for (size_t i = 1; i < numX - 1; ++i) {
            sol.grid(n + 1, i) = sol.grid(n, i) + r * (sol.grid(n, i + 1) - 2.0 * sol.grid(n, i) + sol.grid(n, i - 1));
        }
    }

    return sol;
}

PDESolution PDESolverLab::solveLaplaceEquation2D(
    double Lx,
    double Ly,
    size_t nx,
    size_t ny,
    double top,
    double bottom,
    double left,
    double right,
    double omega,
    double tol,
    int maxIter
) {
    PDESolution sol;
    sol.pdeType = "Elliptic: 2D Laplace Equation (Steady State)";
    sol.method = (omega > 1.0) ? "Gauss-Seidel with Successive Over-Relaxation (SOR)" : "Standard Gauss-Seidel";

    sol.grid = core::Matrix::zeros(ny, nx);
    sol.xCoords = core::Vector::linspace(0.0, Lx, nx);
    sol.yCoordsOrTime = core::Vector::linspace(0.0, Ly, ny);

    // Set boundary conditions
    for (size_t j = 0; j < nx; ++j) {
        sol.grid(0, j) = bottom;        // y = 0
        sol.grid(ny - 1, j) = top;      // y = Ly
    }
    for (size_t i = 0; i < ny; ++i) {
        sol.grid(i, 0) = left;          // x = 0
        sol.grid(i, nx - 1) = right;    // x = Lx
    }

    // Initialize interior to average of boundary values
    double avgBound = 0.25 * (top + bottom + left + right);
    for (size_t i = 1; i < ny - 1; ++i) {
        for (size_t j = 1; j < nx - 1; ++j) {
            sol.grid(i, j) = avgBound;
        }
    }

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        double maxDiff = 0.0;

        for (size_t i = 1; i < ny - 1; ++i) {
            for (size_t j = 1; j < nx - 1; ++j) {
                double uOld = sol.grid(i, j);
                double uGaussSeidel = 0.25 * (sol.grid(i + 1, j) + sol.grid(i - 1, j) +
                                              sol.grid(i, j + 1) + sol.grid(i, j - 1));
                double uNew = (1.0 - omega) * uOld + omega * uGaussSeidel;
                sol.grid(i, j) = uNew;

                maxDiff = std::max(maxDiff, std::abs(uNew - uOld));
            }
        }

        if (maxDiff < tol) {
            break;
        }
    }

    sol.iterations = iter;
    return sol;
}

PDESolution PDESolverLab::solveWaveEquation1D(
    double c,
    double L,
    double tMax,
    size_t numX,
    size_t numT,
    core::RealFunction initialDisplacement,
    core::RealFunction initialVelocity,
    double boundaryLeft,
    double boundaryRight
) {
    PDESolution sol;
    sol.pdeType = "Hyperbolic: 1D Wave Equation";
    sol.method = "Explicit Central Difference Scheme";

    double dx = L / static_cast<double>(numX - 1);
    double dt = tMax / static_cast<double>(numT - 1);
    double C = c * dt / dx; // Courant number

    sol.parameterR = C;
    sol.isStable = (C <= 1.0);

    sol.grid = core::Matrix::zeros(numT, numX);
    sol.xCoords = core::Vector::linspace(0.0, L, numX);
    sol.yCoordsOrTime = core::Vector::linspace(0.0, tMax, numT);

    // Step 0: Initial displacement u(x, 0)
    for (size_t i = 0; i < numX; ++i) {
        sol.grid(0, i) = initialDisplacement(sol.xCoords[i]);
    }
    sol.grid(0, 0) = boundaryLeft;
    sol.grid(0, numX - 1) = boundaryRight;

    // Step 1: Using initial velocity du/dt(x, 0) = g(x)
    // u_i^1 = u_i^0 + dt * g_i + 0.5 * C^2 * (u_{i+1}^0 - 2 u_i^0 + u_{i-1}^0)
    sol.grid(1, 0) = boundaryLeft;
    sol.grid(1, numX - 1) = boundaryRight;
    for (size_t i = 1; i < numX - 1; ++i) {
        double g_i = initialVelocity ? initialVelocity(sol.xCoords[i]) : 0.0;
        sol.grid(1, i) = sol.grid(0, i) + dt * g_i + 0.5 * C * C * (sol.grid(0, i + 1) - 2.0 * sol.grid(0, i) + sol.grid(0, i - 1));
    }

    // Subsequent steps n = 1 to numT - 2:
    // u_i^{n+1} = 2*(1 - C^2)*u_i^n + C^2*(u_{i+1}^n + u_{i-1}^n) - u_i^{n-1}
    for (size_t n = 1; n < numT - 1; ++n) {
        sol.grid(n + 1, 0) = boundaryLeft;
        sol.grid(n + 1, numX - 1) = boundaryRight;

        for (size_t i = 1; i < numX - 1; ++i) {
            sol.grid(n + 1, i) = 2.0 * (1.0 - C * C) * sol.grid(n, i) +
                                 C * C * (sol.grid(n, i + 1) + sol.grid(n, i - 1)) -
                                 sol.grid(n - 1, i);
        }
    }

    return sol;
}

} // namespace ees::unit4

