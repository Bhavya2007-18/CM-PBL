#pragma once

#include <string>
#include <vector>

namespace ees::advisor {

enum class ProblemCategory {
    RootFinding,
    Optimization1D,
    OptimizationMultivariate,
    Interpolation,
    NumericalIntegration,
    LinearSystem,
    Eigenvalues,
    ODEInitialValue,
    PDE
};

struct ProblemProfile {
    ProblemCategory category;
    bool hasBracketingInterval{false};
    bool hasAnalyticDerivative{false};
    bool hasSecondDerivative{false};
    bool isSmoothContinuous{true};
    bool isEquispacedData{true};
    bool isMatrixSymmetricPositiveDefinite{false};
    bool isMatrixSparse{false};
    bool isStiffODE{false};
    int problemDimension{1};
};

struct Recommendation {
    std::string primaryMethod;
    std::string alternateMethod;
    std::vector<std::string> detectedCharacteristics;
    std::string mathematicalJustification;
    std::string cautionaryNotes;
    std::string toString() const;
};

class MethodAdvisor {
public:
    static Recommendation advise(const ProblemProfile& profile);
};

} // namespace ees::advisor

