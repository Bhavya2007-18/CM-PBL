#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include "ees/core/NumericalResult.hpp"
#include <vector>
#include <string>

namespace ees::unit3 {

struct EliminationStep {
    int stepNumber{0};
    std::string operation; // e.g. "Swap R1 <-> R2 (Pivot: 4.5)", "R2 <- R2 - (0.5)*R1"
    core::Matrix currentMatrix;
    core::Vector currentRHS;
};

struct LinearSystemResult {
    bool success{false};
    core::Vector solution;
    core::Matrix finalAugmentedMatrix;
    double determinant{1.0};
    int rowSwaps{0};
    size_t flopCount{0};
    std::string method;
    std::vector<EliminationStep> stepHistory;
    std::string message;
    std::string summary() const;
};

struct FactorizationResult {
    bool success{false};
    core::Matrix L;
    core::Matrix U;
    double reconstructionError{0.0}; // ||A - L*U||_F
    std::string method;
    std::string message;
    std::string summary() const;
};

struct EigenResult {
    bool success{false};
    double eigenvalue{0.0};
    core::Vector eigenvector;
    int iterations{0};
    double finalError{0.0};
    std::vector<core::IterationRecord> history;
    std::string summary() const;
};

class LinearAlgebraLab {
public:
    // Gauss Elimination with Partial Pivoting
    // Computational complexity: 2/3 * n^3 FLOPs
    static LinearSystemResult gaussElimination(
        core::Matrix A,
        core::Vector b,
        bool enablePivoting = true
    );

    // Gauss-Jordan Elimination (computes RREF and solution directly)
    // Computational complexity: n^3 FLOPs
    static LinearSystemResult gaussJordan(
        core::Matrix A,
        core::Vector b
    );

    // Matrix Inverse via Gauss-Jordan [A | I] -> [I | A^-1]
    static core::Matrix inverse(
        core::Matrix A
    );

    // Pivoting Demonstration: Solves ill-conditioned system with and without pivoting
    // Demonstrates numerical catastrophic instability when pivoting is omitted
    static std::string runPivotingDemonstration();

    // Doolittle LU Factorization: A = L * U, where L has 1s on diagonal
    static FactorizationResult doolittleLU(
        const core::Matrix& A
    );

    // Crout LU Factorization: A = L * U, where U has 1s on diagonal
    static FactorizationResult croutLU(
        const core::Matrix& A
    );

    // Cholesky Factorization: A = L * L^T for symmetric positive-definite matrix
    static FactorizationResult cholesky(
        const core::Matrix& A
    );

    // Power Method: Computes dominant eigenvalue and eigenvector
    // Rayleigh quotient iteration: lambda_k = (x_k^T * A * x_k) / (x_k^T * x_k)
    static EigenResult powerMethod(
        const core::Matrix& A,
        core::Vector initialGuess = core::Vector(),
        double tol = 1e-7,
        int maxIter = 100
    );
};

} // namespace ees::unit3

