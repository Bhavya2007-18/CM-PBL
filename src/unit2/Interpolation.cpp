#include "ees/unit2/Interpolation.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace ees::unit2 {

std::string DifferenceTableResult::toString() const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);
    ss << "=== Finite Difference Table ===\n";
    ss << "Number of Points: " << x.size() << "\n";
    if (isEquispaced) {
        ss << "Equispaced interval h = " << h << "\n";
    } else {
        ss << "Unequally spaced data points (Divided Differences applied)\n";
    }

    size_t n = x.size();
    ss << std::setw(10) << "x" << " | " << std::setw(12) << "y";
    for (size_t col = 1; col < n; ++col) {
        ss << " | " << std::setw(12) << ("D^" + std::to_string(col) + " y");
    }
    ss << "\n" << std::string(14 + n * 15, '-') << "\n";

    for (size_t r = 0; r < n; ++r) {
        ss << std::setw(10) << x[r] << " | ";
        for (size_t c = 0; c < n - r; ++c) {
            ss << std::setw(12) << diffTable(r, c) << "   ";
        }
        ss << "\n";
    }
    return ss.str();
}

std::string InterpolationResult::summary() const {
    std::ostringstream ss;
    ss << "=== Interpolation Result (" << method << ") ===\n";
    ss << "Evaluation Point x:  " << std::fixed << std::setprecision(6) << evaluatedX << "\n";
    ss << "Interpolated Value y:" << std::setprecision(8) << interpolatedY << "\n";
    if (estimatedError > 0.0) {
        ss << "Estimated Error:     " << std::scientific << estimatedError << "\n";
    }
    if (!polynomialString.empty()) {
        ss << "Polynomial Form:     " << polynomialString << "\n";
    }
    return ss.str();
}

DifferenceTableResult InterpolationLab::computeForwardDifferenceTable(
    const core::Vector& x,
    const core::Vector& y
) {
    size_t n = x.size();
    if (n != y.size() || n == 0) {
        throw std::invalid_argument("Vectors x and y must have equal positive size");
    }

    DifferenceTableResult res;
    res.x = x;
    res.y = y;
    res.diffTable = core::Matrix::zeros(n, n);

    // Column 0 is y
    for (size_t i = 0; i < n; ++i) {
        res.diffTable(i, 0) = y[i];
    }

    // Check equispaced
    res.isEquispaced = true;
    if (n > 1) {
        res.h = x[1] - x[0];
        for (size_t i = 2; i < n; ++i) {
            if (std::abs((x[i] - x[i - 1]) - res.h) > 1e-7) {
                res.isEquispaced = false;
                break;
            }
        }
    }

    // Compute forward differences: Delta^(col) y_i = Delta^(col-1) y_{i+1} - Delta^(col-1) y_i
    for (size_t col = 1; col < n; ++col) {
        for (size_t row = 0; row < n - col; ++row) {
            res.diffTable(row, col) = res.diffTable(row + 1, col - 1) - res.diffTable(row, col - 1);
        }
    }

    return res;
}

DifferenceTableResult InterpolationLab::computeBackwardDifferenceTable(
    const core::Vector& x,
    const core::Vector& y
) {
    // Backward difference table has the same entries, but referenced from bottom (yn)
    return computeForwardDifferenceTable(x, y);
}

DifferenceTableResult InterpolationLab::computeDividedDifferenceTable(
    const core::Vector& x,
    const core::Vector& y
) {
    size_t n = x.size();
    if (n != y.size() || n == 0) {
        throw std::invalid_argument("Vectors x and y must have equal positive size");
    }

    DifferenceTableResult res;
    res.x = x;
    res.y = y;
    res.isEquispaced = false;
    res.diffTable = core::Matrix::zeros(n, n);

    for (size_t i = 0; i < n; ++i) {
        res.diffTable(i, 0) = y[i];
    }

    // f[x_i, ..., x_{i+k}] = (f[x_{i+1}, ..., x_{i+k}] - f[x_i, ..., x_{i+k-1}]) / (x_{i+k} - x_i)
    for (size_t col = 1; col < n; ++col) {
        for (size_t row = 0; row < n - col; ++row) {
            double denom = x[row + col] - x[row];
            if (std::abs(denom) < 1e-15) {
                throw std::domain_error("Duplicate x values detected in divided difference table");
            }
            res.diffTable(row, col) = (res.diffTable(row + 1, col - 1) - res.diffTable(row, col - 1)) / denom;
        }
    }

    return res;
}

InterpolationResult InterpolationLab::gregoryNewtonForward(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    auto table = computeForwardDifferenceTable(x, y);
    if (!table.isEquispaced) {
        throw std::invalid_argument("Gregory-Newton Forward interpolation requires equispaced data points");
    }

    InterpolationResult res;
    res.method = "Gregory-Newton Forward";
    res.evaluatedX = xVal;

    size_t n = x.size();
    double h = table.h;
    double u = (xVal - x[0]) / h;

    double sum = table.diffTable(0, 0);
    double uTerm = 1.0;
    double factorial = 1.0;

    std::ostringstream poly;
    poly << "P(x) = " << table.diffTable(0, 0);

    for (size_t k = 1; k < n; ++k) {
        uTerm *= (u - static_cast<double>(k - 1));
        factorial *= static_cast<double>(k);
        double term = (uTerm / factorial) * table.diffTable(0, k);
        sum += term;

        poly << " + (" << table.diffTable(0, k) << "/" << k << "!)";
        for (size_t j = 0; j < k; ++j) {
            poly << "*(u - " << j << ")";
        }
    }

    res.interpolatedY = sum;
    res.polynomialString = poly.str() + " [where u = (x - " + std::to_string(x[0]) + ")/" + std::to_string(h) + "]";
    return res;
}

InterpolationResult InterpolationLab::gregoryNewtonBackward(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    auto table = computeBackwardDifferenceTable(x, y);
    if (!table.isEquispaced) {
        throw std::invalid_argument("Gregory-Newton Backward interpolation requires equispaced data points");
    }

    InterpolationResult res;
    res.method = "Gregory-Newton Backward";
    res.evaluatedX = xVal;

    size_t n = x.size();
    double h = table.h;
    size_t lastIdx = n - 1;
    double v = (xVal - x[lastIdx]) / h;

    double sum = table.diffTable(lastIdx, 0);
    double vTerm = 1.0;
    double factorial = 1.0;

    std::ostringstream poly;
    poly << "P(x) = " << table.diffTable(lastIdx, 0);

    for (size_t k = 1; k < n; ++k) {
        vTerm *= (v + static_cast<double>(k - 1));
        factorial *= static_cast<double>(k);
        // In the forward table, the k-th backward difference of y_n is at row (n - 1 - k), column k
        double deltaVal = table.diffTable(lastIdx - k, k);
        double term = (vTerm / factorial) * deltaVal;
        sum += term;

        poly << " + (" << deltaVal << "/" << k << "!)";
        for (size_t j = 0; j < k; ++j) {
            poly << "*(v + " << j << ")";
        }
    }

    res.interpolatedY = sum;
    res.polynomialString = poly.str() + " [where v = (x - " + std::to_string(x[lastIdx]) + ")/" + std::to_string(h) + "]";
    return res;
}

InterpolationResult InterpolationLab::lagrange(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    size_t n = x.size();
    if (n != y.size() || n == 0) {
        throw std::invalid_argument("Vectors x and y must have equal positive size");
    }

    InterpolationResult res;
    res.method = "Lagrange Polynomial";
    res.evaluatedX = xVal;

    double total = 0.0;
    std::ostringstream poly;
    poly << "P(x) = ";

    for (size_t i = 0; i < n; ++i) {
        double Li = 1.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                double denom = x[i] - x[j];
                if (std::abs(denom) < 1e-15) {
                    throw std::domain_error("Duplicate x values detected in Lagrange interpolation");
                }
                Li *= (xVal - x[j]) / denom;
            }
        }
        total += y[i] * Li;

        if (i > 0) poly << " + ";
        poly << y[i] << "*L_" << i << "(x)";
    }

    res.interpolatedY = total;
    res.polynomialString = poly.str();
    return res;
}

InterpolationResult InterpolationLab::newtonDividedDifference(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    auto table = computeDividedDifferenceTable(x, y);

    InterpolationResult res;
    res.method = "Newton Divided Difference";
    res.evaluatedX = xVal;

    size_t n = x.size();
    double sum = table.diffTable(0, 0);
    double xProd = 1.0;

    std::ostringstream poly;
    poly << "P(x) = " << table.diffTable(0, 0);

    for (size_t k = 1; k < n; ++k) {
        xProd *= (xVal - x[k - 1]);
        double coeff = table.diffTable(0, k);
        sum += coeff * xProd;

        poly << " + (" << coeff << ")";
        for (size_t j = 0; j < k; ++j) {
            poly << "*(x - " << x[j] << ")";
        }
    }

    res.interpolatedY = sum;
    res.polynomialString = poly.str();
    return res;
}

std::string InterpolationLab::compareInterpolation(
    const core::Vector& x,
    const core::Vector& y,
    double xVal
) {
    auto gnF = gregoryNewtonForward(x, y, xVal);
    auto gnB = gregoryNewtonBackward(x, y, xVal);
    auto lag = lagrange(x, y, xVal);
    auto ndd = newtonDividedDifference(x, y, xVal);

    std::ostringstream ss;
    ss << "========================================================================================\n";
    ss << "                          INTERPOLATION METHOD COMPARISON                               \n";
    ss << "========================================================================================\n";
    ss << "Evaluation Point: x = " << xVal << "\n";
    ss << std::left << std::setw(30) << "Method"
       << std::setw(20) << "Interpolated y"
       << std::setw(35) << "Domain Suitability" << "\n";
    ss << std::string(85, '-') << "\n";

    ss << std::left << std::setw(30) << "Gregory-Newton Forward"
       << std::setprecision(8) << std::fixed << std::setw(20) << gnF.interpolatedY
       << std::setw(35) << "Near beginning of equispaced table" << "\n";

    ss << std::left << std::setw(30) << "Gregory-Newton Backward"
       << std::setprecision(8) << std::fixed << std::setw(20) << gnB.interpolatedY
       << std::setw(35) << "Near end of equispaced table" << "\n";

    ss << std::left << std::setw(30) << "Lagrange Polynomial"
       << std::setprecision(8) << std::fixed << std::setw(20) << lag.interpolatedY
       << std::setw(35) << "Arbitrary spacing (recomputation on new pts)" << "\n";

    ss << std::left << std::setw(30) << "Newton Divided Difference"
       << std::setprecision(8) << std::fixed << std::setw(20) << ndd.interpolatedY
       << std::setw(35) << "Arbitrary spacing (incremental degree addition)" << "\n";
    ss << std::string(85, '-') << "\n";

    return ss.str();
}

} // namespace ees::unit2

