#include "ees/core/ErrorAnalysis.hpp"
#include <cmath>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <bitset>

namespace ees::core {

double ErrorAnalysis::absoluteError(double trueVal, double approxVal) {
    return std::abs(trueVal - approxVal);
}

double ErrorAnalysis::relativeError(double trueVal, double approxVal) {
    double absTrue = std::abs(trueVal);
    if (absTrue < 1e-30) {
        return std::abs(approxVal);
    }
    return std::abs(trueVal - approxVal) / absTrue;
}

double ErrorAnalysis::percentageError(double trueVal, double approxVal) {
    return relativeError(trueVal, approxVal) * 100.0;
}

double ErrorAnalysis::approximateRelativeError(double currentApprox, double previousApprox) {
    double absCurr = std::abs(currentApprox);
    if (absCurr < 1e-30) {
        return std::abs(currentApprox - previousApprox);
    }
    return (std::abs(currentApprox - previousApprox) / absCurr) * 100.0;
}

int ErrorAnalysis::significantDigits(double relError) {
    if (relError <= 0.0 || !std::isfinite(relError)) {
        return 16;
    }
    double m = -std::floor(std::log10(2.0 * relError));
    if (m < 0) return 0;
    if (m > 16) return 16;
    return static_cast<int>(m);
}

FloatInspection ErrorAnalysis::inspectDouble(double value) {
    FloatInspection info;
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(double));
    info.rawBits = bits;

    info.sign = static_cast<int>((bits >> 63) & 1);
    info.exponent = static_cast<int>((bits >> 52) & 0x7FF);
    info.mantissa = bits & 0xFFFFFFFFFFFFFULL;

    std::string bStr = std::bitset<64>(bits).to_string();
    info.binaryString = bStr.substr(0, 1) + " " + bStr.substr(1, 11) + " " + bStr.substr(12);

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(12) << value;
    info.scientificForm = ss.str();

    std::ostringstream expl;
    expl << "IEEE-754 64-bit Double Precision:\n";
    expl << "  Sign bit:     " << info.sign << " (" << (info.sign ? "Negative" : "Positive") << ")\n";
    expl << "  Exponent:     " << info.exponent << " (Unbiased: " << (info.exponent - 1023) << ")\n";
    expl << "  Mantissa:     0x" << std::hex << info.mantissa << std::dec << " (" << info.mantissa << ")\n";
    info.breakdownExplanation = expl.str();

    return info;
}

FloatInspection ErrorAnalysis::inspectFloat(float value) {
    FloatInspection info;
    uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(float));
    info.rawBits = bits;

    info.sign = static_cast<int>((bits >> 31) & 1);
    info.exponent = static_cast<int>((bits >> 23) & 0xFF);
    info.mantissa = bits & 0x7FFFFFULL;

    std::string bStr = std::bitset<32>(bits).to_string();
    info.binaryString = bStr.substr(0, 1) + " " + bStr.substr(1, 8) + " " + bStr.substr(9);

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6) << value;
    info.scientificForm = ss.str();

    std::ostringstream expl;
    expl << "IEEE-754 32-bit Single Precision:\n";
    expl << "  Sign bit:     " << info.sign << " (" << (info.sign ? "Negative" : "Positive") << ")\n";
    expl << "  Exponent:     " << info.exponent << " (Unbiased: " << (info.exponent - 127) << ")\n";
    expl << "  Mantissa:     0x" << std::hex << info.mantissa << std::dec << " (" << info.mantissa << ")\n";
    info.breakdownExplanation = expl.str();

    return info;
}

double ErrorAnalysis::cancellationDirect(double x) {
    return std::sqrt(x + 1.0) - std::sqrt(x);
}

double ErrorAnalysis::cancellationStable(double x) {
    return 1.0 / (std::sqrt(x + 1.0) + std::sqrt(x));
}

double ErrorAnalysis::cancellationTrigDirect(double x) {
    return (1.0 - std::cos(x)) / (x * x);
}

double ErrorAnalysis::cancellationTrigStable(double x) {
    double halfX = x * 0.5;
    if (std::abs(halfX) < 1e-15) return 0.5;
    double s = std::sin(halfX) / halfX;
    return 0.5 * s * s;
}

} // namespace ees::core
