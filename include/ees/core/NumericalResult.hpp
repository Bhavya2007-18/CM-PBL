#pragma once

#include "Types.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include <vector>
#include <map>
#include <string>

namespace ees::core {

struct IterationRecord {
    int iteration{0};
    double x{0.0};              // Key scalar variable (e.g. current root/estimate)
    Vector xVector;             // Key vector variable (multivariate optimization / power method)
    double fx{0.0};             // Function evaluation f(x)
    double absoluteError{0.0};  // |x_curr - x_prev| or |x_true - x_curr|
    double relativeError{0.0};  // Relative error
    std::map<std::string, double> aux; // Method-specific data (e.g. "a", "b", "f(a)", "f(b)", "f'(x)", "lambda")
    std::string note;           // Row operation note or step explanation
};

struct NumericalResult {
    bool success{false};
    double solution{0.0};
    Vector solutionVector;
    Matrix solutionMatrix;
    int iterations{0};
    double finalError{0.0};
    double finalRelativeError{0.0};
    double executionTimeMs{0.0};
    ConvergenceStatus status{ConvergenceStatus::InvalidInput};
    ConvergenceRate convergenceRate{ConvergenceRate::Unknown};
    std::string convergenceRateExplanation;
    std::string message;
    std::vector<IterationRecord> history;

    // Helper utilities
    void addRecord(const IterationRecord& record) {
        history.push_back(record);
        iterations = static_cast<int>(history.size());
    }

    std::string summary() const;
    std::string tableString(const std::vector<std::string>& auxKeys = {}) const;
};

} // namespace ees::core
