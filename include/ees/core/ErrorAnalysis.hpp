#pragma once

#include "Types.hpp"
#include <cstdint>
#include <string>

namespace ees::core {

struct FloatInspection {
    uint64_t rawBits{0};
    int sign{0};
    int exponent{0};
    uint64_t mantissa{0};
    std::string binaryString;
    std::string scientificForm;
    std::string breakdownExplanation;
};

class ErrorAnalysis {
public:
    // Core error metrics
    static double absoluteError(double trueVal, double approxVal);
    static double relativeError(double trueVal, double approxVal);
    static double percentageError(double trueVal, double approxVal);
    static double approximateRelativeError(double currentApprox, double previousApprox);
    static int significantDigits(double relError);

    // IEEE-754 bit-level representation inspection
    static FloatInspection inspectDouble(double value);
    static FloatInspection inspectFloat(float value);

    // Loss of significance / Catastrophic cancellation demonstration
    // E.g., f(x) = sqrt(x + 1) - sqrt(x) when x is large
    // Direct computation suffers catastrophic cancellation: sqrt(x+1) - sqrt(x)
    // Mathematically equivalent stable form: 1 / (sqrt(x+1) + sqrt(x))
    static double cancellationDirect(double x);
    static double cancellationStable(double x);

    // E.g., (1 - cos(x)) / (x^2) for small x
    static double cancellationTrigDirect(double x);
    static double cancellationTrigStable(double x);
};

} // namespace ees::core
