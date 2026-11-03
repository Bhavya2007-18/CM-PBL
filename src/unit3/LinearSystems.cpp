#include "ees/unit3/LinearSystems.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace ees::unit3 {

std::string LinearSystemResult::summary() const {
    std::ostringstream ss;
    ss << "=== Linear System Solution (" << method << ") ===\n";
    ss << "Status: " << (success ? "SOLVED" : "FAILED") << " - " << message << "\n";
    if (success) {
        ss << "Solution Vector x:\n" << solution.toString(8) << "\n";
        ss << "Determinant of A: " << std::setprecision(6) << determinant << "\n";
        ss << "Row Swaps:        " << rowSwaps << "\n";
        ss << "Estimated FLOPs:  " << flopCount << "\n";
    }
    return ss.str();
}

std::string FactorizationResult::summary() const {
    std::ostringstream ss;
    ss << "=== Matrix Factorization (" << method << ") ===\n";
    ss << "Status: " << (success ? "SUCCESS" : "FAILED") << "\n";
    if (success) {
        ss << "Lower Matrix L:\n" << L.toString(6) << "\n\n";
        ss << "Upper Matrix U:\n" << U.toString(6) << "\n\n";
        ss << "Reconstruction Error ||A - LU||_F: " << std::scientific << reconstructionError << "\n";
    } else {
        ss << "Message: " << message << "\n";
    }
    return ss.str();
}

std::string EigenResult::summary() const {
    std::ostringstream ss;
    ss << "=== Power Method Eigenvalue Analysis ===\n";
    ss << "Status: " << (success ? "CONVERGED" : "FAILED") << "\n";
    if (success) {
        ss << "Dominant Eigenvalue (lambda): " << std::setprecision(8) << eigenvalue << "\n";
        ss << "Dominant Eigenvector:\n" << eigenvector.toString(6) << "\n";
        ss << "Iterations: " << iterations << "\n";
        ss << "Final Error: " << std::scientific << finalError << "\n";
    }
    return ss.str();
}

LinearSystemResult LinearAlgebraLab::gaussElimination(
    core::Matrix A,
    core::Vector b,
    bool enablePivoting
) {
    LinearSystemResult res;
    res.method = enablePivoting ? "Gauss Elimination with Partial Pivoting" : "Naive Gauss Elimination";

    const size_t n = A.rows();
    if (!A.isSquare() || b.size() != n) {
        res.success = false;
        res.message = "System dimension mismatch: A must be n x n and b must have length n";
        return res;
    }

    size_t flops = 0;
    int sign = 1;
    res.stepHistory.push_back({0, "Initial System", A, b});

    // Forward Elimination
    for (size_t k = 0; k < n - 1; ++k) {
        if (enablePivoting) {
            // Find pivot
            size_t maxRow = k;
            double maxVal = std::abs(A(k, k));
            for (size_t i = k + 1; i < n; ++i) {
                if (std::abs(A(i, k)) > maxVal) {
                    maxVal = std::abs(A(i, k));
                    maxRow = i;
                }
            }

            if (maxVal < 1e-14) {
                res.success = false;
                res.message = "Matrix is singular or near-singular: zero pivot encountered at column " + std::to_string(k);
                return res;
            }

            if (maxRow != k) {
                A.swapRows(k, maxRow);
                std::swap(b[k], b[maxRow]);
                sign = -sign;
                res.rowSwaps++;

                std::string op = "Swap R" + std::to_string(k + 1) + " <-> R" + std::to_string(maxRow + 1) +
                                 " (Pivot: " + std::to_string(A(k, k)) + ")";
                res.stepHistory.push_back({static_cast<int>(res.stepHistory.size()), op, A, b});
            }
        } else {
            if (std::abs(A(k, k)) < 1e-14) {
                res.success = false;
                res.message = "Zero pivot encountered in naive elimination at row " + std::to_string(k);
                return res;
            }
        }

        // Eliminate entries below pivot
        for (size_t i = k + 1; i < n; ++i) {
            double factor = A(i, k) / A(k, k);
            flops += 1; // division

            for (size_t j = k; j < n; ++j) {
                A(i, j) -= factor * A(k, j);
                flops += 2; // mul + sub
            }
            b[i] -= factor * b[k];
            flops += 2; // mul + sub

            std::ostringstream ss;
            ss << "R" << (i + 1) << " <- R" << (i + 1) << " - (" << factor << ")*R" << (k + 1);
            res.stepHistory.push_back({static_cast<int>(res.stepHistory.size()), ss.str(), A, b});
        }
    }

    if (std::abs(A(n - 1, n - 1)) < 1e-14) {
        res.success = false;
        res.message = "Singular matrix: zero element in last diagonal entry";
        return res;
    }

    // Compute determinant
    double det = sign;
    for (size_t i = 0; i < n; ++i) det *= A(i, i);
    res.determinant = det;

    // Back Substitution
    core::Vector x(n, 0.0);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        double sum = b[i];
        for (size_t j = i + 1; j < n; ++j) {
            sum -= A(i, j) * x[j];
            flops += 2;
        }
        x[i] = sum / A(i, i);
        flops += 1;
    }

    res.success = true;
    res.solution = x;
    res.finalAugmentedMatrix = A.augment(b);
    res.flopCount = flops;
    res.message = "System solved successfully.";
    return res;
}

LinearSystemResult LinearAlgebraLab::gaussJordan(
    core::Matrix A,
    core::Vector b
) {
    LinearSystemResult res;
    res.method = "Gauss-Jordan Elimination (RREF)";

    const size_t n = A.rows();
    if (!A.isSquare() || b.size() != n) {
        res.success = false;
        res.message = "System dimension mismatch";
        return res;
    }

    size_t flops = 0;
    res.stepHistory.push_back({0, "Initial System", A, b});

    for (size_t k = 0; k < n; ++k) {
        // Partial pivoting
        size_t maxRow = k;
        double maxVal = std::abs(A(k, k));
        for (size_t i = k + 1; i < n; ++i) {
            if (std::abs(A(i, k)) > maxVal) {
                maxVal = std::abs(A(i, k));
                maxRow = i;
            }
        }

        if (maxVal < 1e-14) {
            res.success = false;
            res.message = "Matrix is singular: no non-zero pivot found in column " + std::to_string(k);
            return res;
        }

        if (maxRow != k) {
            A.swapRows(k, maxRow);
            std::swap(b[k], b[maxRow]);
            res.rowSwaps++;
            res.stepHistory.push_back({static_cast<int>(res.stepHistory.size()),
                "Swap R" + std::to_string(k + 1) + " <-> R" + std::to_string(maxRow + 1), A, b});
        }

        // Scale pivot row so A(k, k) == 1
        double pivot = A(k, k);
        for (size_t j = k; j < n; ++j) {
            A(k, j) /= pivot;
            flops += 1;
        }
        b[k] /= pivot;
        flops += 1;

        res.stepHistory.push_back({static_cast<int>(res.stepHistory.size()),
            "Scale R" + std::to_string(k + 1) + " by 1 / (" + std::to_string(pivot) + ")", A, b});

        // Eliminate all other rows (both above and below pivot)
        for (size_t i = 0; i < n; ++i) {
            if (i != k) {
                double factor = A(i, k);
                for (size_t j = k; j < n; ++j) {
                    A(i, j) -= factor * A(k, j);
                    flops += 2;
                }
                b[i] -= factor * b[k];
                flops += 2;
            }
        }
        res.stepHistory.push_back({static_cast<int>(res.stepHistory.size()),
            "Eliminate column " + std::to_string(k + 1) + " across all other rows", A, b});
    }

    res.success = true;
    res.solution = b; // In RREF, the RHS contains the solution directly!
    res.finalAugmentedMatrix = A.augment(b);
    res.flopCount = flops;
    res.message = "System transformed into Reduced Row Echelon Form (RREF).";
    return res;
}

core::Matrix LinearAlgebraLab::inverse(core::Matrix A) {
    if (!A.isSquare()) throw std::invalid_argument("Matrix must be square to invert");
    const size_t n = A.rows();

    core::Matrix aug = A.augment(core::Matrix::identity(n));

    for (size_t k = 0; k < n; ++k) {
        size_t maxRow = k;
        double maxVal = std::abs(aug(k, k));
        for (size_t i = k + 1; i < n; ++i) {
            if (std::abs(aug(i, k)) > maxVal) {
                maxVal = std::abs(aug(i, k));
                maxRow = i;
            }
        }
        if (maxVal < 1e-14) throw std::runtime_error("Matrix is singular: cannot compute inverse");

        if (maxRow != k) aug.swapRows(k, maxRow);

        double pivot = aug(k, k);
        for (size_t j = k; j < 2 * n; ++j) {
            aug(k, j) /= pivot;
        }

        for (size_t i = 0; i < n; ++i) {
            if (i != k) {
                double factor = aug(i, k);
                for (size_t j = k; j < 2 * n; ++j) {
                    aug(i, j) -= factor * aug(k, j);
                }
            }
        }
    }

    return aug.subMatrix(0, n, n, n);
}

std::string LinearAlgebraLab::runPivotingDemonstration() {
    std::ostringstream ss;
    ss << "========================================================================\n";
    ss << "                  PIVOTING & NUMERICAL STABILITY LAB                    \n";
    ss << "========================================================================\n";
    ss << "Classic Wilkinson/Forsythe Ill-Conditioned System:\n";
    ss << "  [ eps   1.0 ] [ x1 ] = [ 1.0 ]\n";
    ss << "  [ 1.0   1.0 ] [ x2 ] = [ 2.0 ]\n";
    ss << "where eps = 1e-12 (a tiny number compared to 1.0).\n\n";

    double eps = 1e-12;
    core::Matrix A{
        {eps, 1.0},
        {1.0, 1.0}
    };
    core::Vector b{1.0, 2.0};

    auto resNaive = gaussElimination(A, b, false);
    auto resPiv   = gaussElimination(A, b, true);

    ss << "1. WITHOUT PIVOTING (Naive Elimination):\n";
    ss << "   Computed x1 = " << std::setprecision(10) << resNaive.solution[0] << "\n";
    ss << "   Computed x2 = " << resNaive.solution[1] << "\n\n";

    ss << "2. WITH PARTIAL PIVOTING:\n";
    ss << "   Computed x1 = " << resPiv.solution[0] << "\n";
    ss << "   Computed x2 = " << resPiv.solution[1] << "\n\n";

    ss << "Exact Analytical Solution: x1 = 1 / (1 - eps) ≈ 1.000000000001, x2 ≈ 0.999999999999.\n";
    ss << "Analysis: Without pivoting, dividing by tiny eps = 1e-12 causes catastrophic\n";
    ss << "swamping of significant digits in standard floating-point representation.\n";
    ss << "Partial pivoting eliminates this instability by always swapping the largest\n";
    ss << "magnitude element into the pivot position.\n";

    return ss.str();
}

FactorizationResult LinearAlgebraLab::doolittleLU(const core::Matrix& A) {
    FactorizationResult res;
    res.method = "Doolittle LU Factorization (Unit Lower Triangular)";

    if (!A.isSquare()) {
        res.success = false;
        res.message = "Matrix must be square for LU factorization.";
        return res;
    }

    const size_t n = A.rows();
    core::Matrix L = core::Matrix::identity(n);
    core::Matrix U = core::Matrix::zeros(n, n);

    for (size_t i = 0; i < n; ++i) {
        // Upper Triangular U
        for (size_t k = i; k < n; ++k) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += L(i, j) * U(j, k);
            }
            U(i, k) = A(i, k) - sum;
        }

        if (std::abs(U(i, i)) < 1e-14 && i + 1 < n) {
            res.success = false;
            res.message = "Zero pivot on U diagonal: Doolittle factorization requires pivoting.";
            return res;
        }

        // Lower Triangular L
        for (size_t k = i + 1; k < n; ++k) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += L(k, j) * U(j, i);
            }
            L(k, i) = (A(k, i) - sum) / U(i, i);
        }
    }

    res.success = true;
    res.L = L;
    res.U = U;
    res.reconstructionError = (A - L * U).normFrobenius();
    return res;
}

FactorizationResult LinearAlgebraLab::croutLU(const core::Matrix& A) {
    FactorizationResult res;
    res.method = "Crout LU Factorization (Unit Upper Triangular)";

    if (!A.isSquare()) {
        res.success = false;
        res.message = "Matrix must be square.";
        return res;
    }

    const size_t n = A.rows();
    core::Matrix L = core::Matrix::zeros(n, n);
    core::Matrix U = core::Matrix::identity(n);

    for (size_t j = 0; j < n; ++j) {
        // Lower Triangular L column j
        for (size_t i = j; i < n; ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < j; ++k) {
                sum += L(i, k) * U(k, j);
            }
            L(i, j) = A(i, j) - sum;
        }

        if (std::abs(L(j, j)) < 1e-14 && j + 1 < n) {
            res.success = false;
            res.message = "Zero diagonal entry in L.";
            return res;
        }

        // Upper Triangular U row j
        for (size_t i = j + 1; i < n; ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < j; ++k) {
                sum += L(j, k) * U(k, i);
            }
            U(j, i) = (A(j, i) - sum) / L(j, j);
        }
    }

    res.success = true;
    res.L = L;
    res.U = U;
    res.reconstructionError = (A - L * U).normFrobenius();
    return res;
}

FactorizationResult LinearAlgebraLab::cholesky(const core::Matrix& A) {
    FactorizationResult res;
    res.method = "Cholesky Factorization (A = L * L^T)";

    if (!A.isSquare()) {
        res.success = false;
        res.message = "Matrix must be square.";
        return res;
    }

    const size_t n = A.rows();
    core::Matrix L = core::Matrix::zeros(n, n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < j; ++k) {
                sum += L(i, k) * L(j, k);
            }

            if (i == j) {
                double val = A(i, i) - sum;
                if (val <= 0.0) {
                    res.success = false;
                    res.message = "Matrix is not positive-definite: negative or zero diagonal pivot at row " + std::to_string(i);
                    return res;
                }
                L(i, i) = std::sqrt(val);
            } else {
                L(i, j) = (A(i, j) - sum) / L(j, j);
            }
        }
    }

    res.success = true;
    res.L = L;
    res.U = L.transpose();
    res.reconstructionError = (A - L * res.U).normFrobenius();
    return res;
}

EigenResult LinearAlgebraLab::powerMethod(
    const core::Matrix& A,
    core::Vector initialGuess,
    double tol,
    int maxIter
) {
    EigenResult res;
    if (!A.isSquare()) {
        res.success = false;
        return res;
    }

    const size_t n = A.rows();
    core::Vector x = initialGuess.empty() ? core::Vector::ones(n) : initialGuess;
    x = x.normalized();

    double prevLambda = 0.0;

    for (int iter = 1; iter <= maxIter; ++iter) {
        core::Vector y = A * x;
        double lambda = x.dot(y); // Rayleigh quotient: x^T A x / (x^T x) since ||x|| = 1

        double err = std::abs(lambda - prevLambda);
        core::IterationRecord rec;
        rec.iteration = iter;
        rec.x = lambda;
        rec.xVector = x;
        rec.absoluteError = err;
        rec.aux["eigenvalue"] = lambda;
        res.history.push_back(rec);

        if (err < tol && iter > 1) {
            res.success = true;
            res.eigenvalue = lambda;
            res.eigenvector = x;
            res.iterations = iter;
            res.finalError = err;
            break;
        }

        if (y.norm2() < 1e-15) {
            res.success = false;
            break;
        }

        x = y.normalized();
        prevLambda = lambda;
    }

    if (!res.success && !res.history.empty()) {
        res.eigenvalue = prevLambda;
        res.eigenvector = x;
        res.iterations = static_cast<int>(res.history.size());
        res.finalError = res.history.back().absoluteError;
    }

    return res;
}

} // namespace ees::unit3

