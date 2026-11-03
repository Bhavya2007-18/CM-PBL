#include "ees/core/NumericalResult.hpp"
#include <sstream>
#include <iomanip>

namespace ees::core {

std::string NumericalResult::summary() const {
    std::ostringstream ss;
    ss << "=== Numerical Computation Result ===\n";
    ss << "Status:           " << statusToString(status) << (success ? " [SUCCESS]" : " [FAILED]") << "\n";
    if (solutionVector.empty() && solutionMatrix.empty()) {
        ss << "Scalar Solution:  " << std::setprecision(10) << solution << "\n";
    } else if (!solutionVector.empty()) {
        ss << "Vector Solution:  " << solutionVector.toString(8) << "\n";
    }
    ss << "Iterations:       " << iterations << "\n";
    ss << "Final Error:      " << std::scientific << std::setprecision(6) << finalError << "\n";
    ss << "Relative Error:   " << std::scientific << std::setprecision(6) << finalRelativeError << "\n";
    ss << "Convergence Rate: " << convergenceRateToString(convergenceRate) << "\n";
    if (!convergenceRateExplanation.empty()) {
        ss << "Rate Analysis:    " << convergenceRateExplanation << "\n";
    }
    ss << "Elapsed Time:     " << std::fixed << std::setprecision(3) << executionTimeMs << " ms\n";
    if (!message.empty()) {
        ss << "Message:          " << message << "\n";
    }
    return ss.str();
}

std::string NumericalResult::tableString(const std::vector<std::string>& auxKeys) const {
    if (history.empty()) return "No iteration records available.\n";

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);

    // Header
    ss << std::setw(6) << "Iter" << " | "
       << std::setw(14) << "x" << " | "
       << std::setw(14) << "f(x)" << " | "
       << std::setw(14) << "Abs Error" << " | "
       << std::setw(14) << "Rel Error";

    for (const auto& key : auxKeys) {
        ss << " | " << std::setw(12) << key;
    }
    ss << "\n";

    // Separator line
    size_t lineLength = 6 + 3 + 14 + 3 + 14 + 3 + 14 + 3 + 14 + auxKeys.size() * 15;
    ss << std::string(lineLength, '-') << "\n";

    // Rows
    for (const auto& rec : history) {
        ss << std::setw(6) << rec.iteration << " | "
           << std::setw(14) << rec.x << " | "
           << std::setw(14) << rec.fx << " | "
           << std::setw(14) << rec.absoluteError << " | "
           << std::setw(14) << rec.relativeError;

        for (const auto& key : auxKeys) {
            auto it = rec.aux.find(key);
            if (it != rec.aux.end()) {
                ss << " | " << std::setw(12) << it->second;
            } else {
                ss << " | " << std::setw(12) << "-";
            }
        }
        if (!rec.note.empty()) {
            ss << "  (" << rec.note << ")";
        }
        ss << "\n";
    }

    return ss.str();
}

} // namespace ees::core
