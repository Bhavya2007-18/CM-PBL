#pragma once

#include <string>
#include <functional>
#include <cstddef>
#include <limits>
#include <cmath>

namespace ees::core {

using Real = double;
using Index = std::size_t;

constexpr double EPSILON = 1e-14;
constexpr double PI_VAL = 3.141592653589793238462643383279502884;
constexpr double E_VAL  = 2.718281828459045235360287471352662497;

enum class ConvergenceStatus {
    Converged,
    MaxIterationsReached,
    Diverged,
    InvalidInput,
    SingularMatrix,
    DivisionByZero,
    LossOfSignificance,
    UserAborted
};

inline std::string statusToString(ConvergenceStatus status) {
    switch (status) {
        case ConvergenceStatus::Converged:             return "Converged";
        case ConvergenceStatus::MaxIterationsReached: return "Max Iterations Reached";
        case ConvergenceStatus::Diverged:             return "Diverged";
        case ConvergenceStatus::InvalidInput:         return "Invalid Input";
        case ConvergenceStatus::SingularMatrix:       return "Singular Matrix";
        case ConvergenceStatus::DivisionByZero:       return "Division by Zero";
        case ConvergenceStatus::LossOfSignificance:   return "Loss of Significance";
        case ConvergenceStatus::UserAborted:          return "User Aborted";
        default:                                      return "Unknown Status";
    }
}

enum class OptimizationGoal {
    Minimize,
    Maximize
};

enum class ConvergenceRate {
    Unknown,
    Sublinear,
    Linear,
    Superlinear,
    Quadratic,
    Cubic
};

inline std::string convergenceRateToString(ConvergenceRate rate) {
    switch (rate) {
        case ConvergenceRate::Sublinear:   return "Sublinear";
        case ConvergenceRate::Linear:      return "Linear (Order 1)";
        case ConvergenceRate::Superlinear: return "Superlinear (Order ~1.618)";
        case ConvergenceRate::Quadratic:   return "Quadratic (Order 2)";
        case ConvergenceRate::Cubic:       return "Cubic (Order 3)";
        default:                           return "Indeterminate / Empirical";
    }
}

// Function pointer/lambda type definitions
using RealFunction = std::function<double(double)>;
using RealFunction2D = std::function<double(double, double)>;

} // namespace ees::core
