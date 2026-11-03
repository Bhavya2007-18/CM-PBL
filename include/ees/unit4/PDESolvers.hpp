#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include <vector>
#include <string>

namespace ees::unit4 {

struct PDESolution {
    core::Matrix grid; // Heat: rows = time steps, cols = spatial points
                       // Laplace: 2D steady-state temperature field Nx x Ny
                       // Wave: rows = time steps, cols = spatial points
    core::Vector xCoords;
    core::Vector yCoordsOrTime;
    std::string pdeType;
    std::string method;
    bool isStable{true};
    double parameterR{0.0}; // CFL stability parameter
    int iterations{0};
    std::string summary() const;
    std::string asciiHeatmap(int displayCols = 40, int displayRows = 15) const;
};

class PDESolverLab {
public:
    // Parabolic PDE: 1D Heat / Diffusion Equation
    // du/dt = alpha * d^2u/dx^2, for x in [0, L], t in [0, tMax]
    // Stability condition: r = alpha * dt / (dx^2) <= 0.5
    static PDESolution solveHeatEquation1D(
        double alpha,
        double L,
        double tMax,
        size_t numX,
        size_t numT,
        core::RealFunction initialCondition, // u(x, 0)
        double boundaryLeft = 0.0,           // u(0, t)
        double boundaryRight = 0.0           // u(L, t)
    );

    // Elliptic PDE: 2D Laplace Equation
    // d^2u/dx^2 + d^2u/dy^2 = 0 on [0, Lx] x [0, Ly]
    // Dirichlet boundary conditions: top, bottom, left, right
    // Solved via Gauss-Seidel with Successive Over-Relaxation (SOR)
    static PDESolution solveLaplaceEquation2D(
        double Lx,
        double Ly,
        size_t nx,
        size_t ny,
        double top,
        double bottom,
        double left,
        double right,
        double omega = 1.25, // SOR relaxation factor (1.0 = Gauss-Seidel)
        double tol = 1e-5,
        int maxIter = 1000
    );

    // Hyperbolic PDE: 1D Wave Equation
    // d^2u/dt^2 = c^2 * d^2u/dx^2, for x in [0, L], t in [0, tMax]
    // Courant-Friedrichs-Lewy (CFL) stability condition: C = c * dt / dx <= 1.0
    static PDESolution solveWaveEquation1D(
        double c,
        double L,
        double tMax,
        size_t numX,
        size_t numT,
        core::RealFunction initialDisplacement, // u(x, 0)
        core::RealFunction initialVelocity,     // du/dt(x, 0)
        double boundaryLeft = 0.0,
        double boundaryRight = 0.0
    );
};

} // namespace ees::unit4

