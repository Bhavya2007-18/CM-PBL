#include "ees/ui/Visualizer.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <fstream>

namespace ees::ui {

std::string Visualizer::plotFunction(
    core::RealFunction f,
    double xMin,
    double xMax,
    int width,
    int height,
    const std::vector<double>& highlightX,
    const std::string& title
) {
    if (width < 20) width = 20;
    if (height < 10) height = 10;
    if (xMax <= xMin) xMax = xMin + 1.0;

    // Sample function across width
    std::vector<double> xVals(width);
    std::vector<double> yVals(width);

    double dx = (xMax - xMin) / (width - 1);
    double yMin = 1e30;
    double yMax = -1e30;

    for (int col = 0; col < width; ++col) {
        double x = xMin + col * dx;
        xVals[col] = x;
        double y = f(x);
        yVals[col] = y;
        if (std::isfinite(y)) {
            yMin = std::min(yMin, y);
            yMax = std::max(yMax, y);
        }
    }

    if (yMax <= yMin) {
        yMax = yMin + 1.0;
    }
    // Add margin
    double yMargin = 0.05 * (yMax - yMin);
    yMin -= yMargin;
    yMax += yMargin;

    // Canvas matrix: height rows, width cols
    std::vector<std::string> canvas(height, std::string(width, ' '));

    // Draw zero-axes if in range
    int zeroRow = -1;
    if (yMin <= 0.0 && yMax >= 0.0) {
        zeroRow = static_cast<int>((yMax / (yMax - yMin)) * (height - 1));
        zeroRow = std::clamp(zeroRow, 0, height - 1);
        for (int c = 0; c < width; ++c) {
            canvas[zeroRow][c] = '-';
        }
    }

    int zeroCol = -1;
    if (xMin <= 0.0 && xMax >= 0.0) {
        zeroCol = static_cast<int>(((-xMin) / (xMax - xMin)) * (width - 1));
        zeroCol = std::clamp(zeroCol, 0, width - 1);
        for (int r = 0; r < height; ++r) {
            canvas[r][zeroCol] = '|';
        }
    }

    if (zeroRow >= 0 && zeroCol >= 0) {
        canvas[zeroRow][zeroCol] = '+';
    }

    // Plot function curve
    for (int c = 0; c < width; ++c) {
        double y = yVals[c];
        if (std::isfinite(y)) {
            int r = static_cast<int>(((yMax - y) / (yMax - yMin)) * (height - 1));
            r = std::clamp(r, 0, height - 1);
            canvas[r][c] = '*';
        }
    }

    // Plot highlight markers (roots, iterations)
    for (double hx : highlightX) {
        if (hx >= xMin && hx <= xMax) {
            int c = static_cast<int>(((hx - xMin) / (xMax - xMin)) * (width - 1));
            double hy = f(hx);
            int r = static_cast<int>(((yMax - hy) / (yMax - yMin)) * (height - 1));
            r = std::clamp(r, 0, height - 1);
            c = std::clamp(c, 0, width - 1);
            canvas[r][c] = 'O'; // Root / Marker
        }
    }

    std::ostringstream ss;
    ss << "=== " << title << " ===\n";
    ss << "  y_max = " << std::fixed << std::setprecision(4) << yMax << "\n";
    for (int r = 0; r < height; ++r) {
        ss << "  |" << canvas[r] << "|\n";
    }
    ss << "  y_min = " << yMin << "\n";
    ss << "  x: [" << xMin << " ... " << xMax << "]\n";
    ss << "  Legend: [*] f(x) curve   [O] Root / Iteration Marker   [+] Origin\n";
    return ss.str();
}

std::string Visualizer::plotConvergence(
    const std::vector<core::IterationRecord>& history,
    int width,
    int height,
    const std::string& title
) {
    if (history.empty()) return "No convergence history to plot.\n";

    size_t numIters = history.size();
    std::vector<double> iters(numIters);
    std::vector<double> logErrors(numIters);

    double minLog = 100.0;
    double maxLog = -100.0;

    for (size_t i = 0; i < numIters; ++i) {
        iters[i] = history[i].iteration;
        double err = history[i].absoluteError;
        double lErr = (err > 1e-18) ? std::log10(err) : -18.0;
        logErrors[i] = lErr;
        minLog = std::min(minLog, lErr);
        maxLog = std::max(maxLog, lErr);
    }

    if (maxLog <= minLog) maxLog = minLog + 1.0;

    std::vector<std::string> canvas(height, std::string(width, ' '));

    for (size_t i = 0; i < numIters; ++i) {
        int c = (numIters > 1) ? static_cast<int>((static_cast<double>(i) / (numIters - 1)) * (width - 1)) : 0;
        int r = static_cast<int>(((maxLog - logErrors[i]) / (maxLog - minLog)) * (height - 1));
        r = std::clamp(r, 0, height - 1);
        c = std::clamp(c, 0, width - 1);
        canvas[r][c] = '#';
    }

    std::ostringstream ss;
    ss << "=== " << title << " ===\n";
    ss << "  log10(Error) [Top = 10^" << std::setprecision(1) << maxLog << "]\n";
    for (int r = 0; r < height; ++r) {
        ss << "  |" << canvas[r] << "|\n";
    }
    ss << "  log10(Error) [Bottom = 10^" << minLog << "]\n";
    ss << "  Iteration: 1 -> " << numIters << "\n";
    return ss.str();
}

std::string Visualizer::plotIntervalShrinkage(
    const std::vector<core::IterationRecord>& history,
    int barWidth
) {
    if (history.empty()) return "No iteration records.\n";

    // Find global [a0, b0] from first record aux
    double globalA = 0.0, globalB = 1.0;
    if (history[0].aux.count("a") && history[0].aux.count("b")) {
        globalA = history[0].aux.at("a");
        globalB = history[0].aux.at("b");
    }

    double globalSpan = globalB - globalA;
    if (std::abs(globalSpan) < 1e-12) globalSpan = 1.0;

    std::ostringstream ss;
    ss << "=== Interval Shrinkage Trajectory ===\n";
    ss << "Iter | Interval [a, b] Bar (" << globalA << " to " << globalB << ") | Width\n";
    ss << std::string(barWidth + 30, '-') << "\n";

    for (const auto& rec : history) {
        if (!rec.aux.count("a") || !rec.aux.count("b")) continue;
        double a = rec.aux.at("a");
        double b = rec.aux.at("b");
        double width = std::abs(b - a);

        int left = static_cast<int>(((a - globalA) / globalSpan) * barWidth);
        int right = static_cast<int>(((b - globalA) / globalSpan) * barWidth);
        left = std::clamp(left, 0, barWidth - 1);
        right = std::clamp(right, 0, barWidth - 1);
        if (right <= left) right = left + 1;

        std::string bar(barWidth, ' ');
        for (int c = left; c <= right && c < barWidth; ++c) {
            bar[c] = '=';
        }
        bar[left] = '[';
        if (right < barWidth) bar[right] = ']';

        ss << std::setw(4) << rec.iteration << " |" << bar << "| "
           << std::scientific << std::setprecision(3) << width << "\n";
    }
    return ss.str();
}

std::string Visualizer::plotMultiSeries(
    const std::vector<PlotSeries>& series,
    int width,
    int height,
    const std::string& title
) {
    if (series.empty()) return "No series provided.\n";

    double xMin = 1e30, xMax = -1e30;
    double yMin = 1e30, yMax = -1e30;

    for (const auto& s : series) {
        for (size_t i = 0; i < s.x.size(); ++i) {
            xMin = std::min(xMin, s.x[i]);
            xMax = std::max(xMax, s.x[i]);
            yMin = std::min(yMin, s.y[i]);
            yMax = std::max(yMax, s.y[i]);
        }
    }

    if (xMax <= xMin) xMax = xMin + 1.0;
    if (yMax <= yMin) yMax = yMin + 1.0;

    std::vector<std::string> canvas(height, std::string(width, ' '));

    for (const auto& s : series) {
        for (size_t i = 0; i < s.x.size(); ++i) {
            int c = static_cast<int>(((s.x[i] - xMin) / (xMax - xMin)) * (width - 1));
            int r = static_cast<int>(((yMax - s.y[i]) / (yMax - yMin)) * (height - 1));
            r = std::clamp(r, 0, height - 1);
            c = std::clamp(c, 0, width - 1);
            canvas[r][c] = s.pointChar;
        }
    }

    std::ostringstream ss;
    ss << "=== " << title << " ===\n";
    ss << "  y_max = " << std::fixed << std::setprecision(4) << yMax << "\n";
    for (int r = 0; r < height; ++r) {
        ss << "  |" << canvas[r] << "|\n";
    }
    ss << "  y_min = " << yMin << "\n";
    ss << "  x: [" << xMin << " ... " << xMax << "]\n";
    ss << "  Series Legend: ";
    for (const auto& s : series) {
        ss << "[" << s.pointChar << "] " << s.label << "   ";
    }
    ss << "\n";
    return ss.str();
}

bool Visualizer::exportToHTML(
    const std::string& filePath,
    const std::string& pageTitle,
    const core::NumericalResult& result,
    const std::string& mathExplanation
) {
    std::ofstream out(filePath);
    if (!out.is_open()) return false;

    out << "<!DOCTYPE html>\n";
    out << "<html lang=\"en\">\n<head>\n";
    out << "<meta charset=\"UTF-8\">\n";
    out << "<title>" << pageTitle << " - EES Dashboard</title>\n";
    out << "<style>\n";
    out << "  body { background-color: #121417; color: #E0E2E6; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; margin: 0; padding: 24px; }\n";
    out << "  h1, h2, h3 { color: #FFFFFF; font-weight: 500; margin-top: 0; }\n";
    out << "  .card { background: #1A1D21; border: 1px solid #2B3037; border-radius: 6px; padding: 20px; margin-bottom: 20px; }\n";
    out << "  .accent { color: #F59E0B; font-weight: 600; }\n";
    out << "  .success { color: #10B981; font-weight: 600; }\n";
    out << "  table { width: 100%; border-collapse: collapse; font-family: 'Consolas', 'Courier New', monospace; font-size: 13px; }\n";
    out << "  th, td { border: 1px solid #2B3037; padding: 8px 12px; text-align: right; }\n";
    out << "  th { background-color: #22262B; color: #9CA3AF; text-transform: uppercase; font-size: 11px; letter-spacing: 0.05em; }\n";
    out << "  tr:nth-child(even) { background-color: #16181B; }\n";
    out << "  pre { background-color: #0E1012; border: 1px solid #22262B; padding: 14px; border-radius: 4px; overflow-x: auto; color: #38BDF8; font-family: 'Consolas', monospace; }\n";
    out << "</style>\n</head>\n<body>\n";

    out << "<h1>EES Scientific Dashboard</h1>\n";
    out << "<h2>" << pageTitle << "</h2>\n";

    // Summary Card
    out << "<div class=\"card\">\n";
    out << "  <h3>Computation Summary</h3>\n";
    out << "  <p>Status: <span class=\"" << (result.success ? "success" : "accent") << "\">"
        << core::statusToString(result.status) << "</span></p>\n";
    out << "  <p>Iterations: <strong>" << result.iterations << "</strong> | Execution Time: <strong>"
        << std::fixed << std::setprecision(3) << result.executionTimeMs << " ms</strong></p>\n";
    out << "  <p>Convergence Rate: <span class=\"accent\">" << core::convergenceRateToString(result.convergenceRate)
        << "</span> (" << result.convergenceRateExplanation << ")</p>\n";
    out << "</div>\n";

    // Mathematical Explanation Card
    if (!mathExplanation.empty()) {
        out << "<div class=\"card\">\n";
        out << "  <h3>Mathematical Intuition & Theory</h3>\n";
        out << "  <pre>" << mathExplanation << "</pre>\n";
        out << "</div>\n";
    }

    // Iteration Table
    out << "<div class=\"card\">\n";
    out << "  <h3>Iteration History Data</h3>\n";
    out << "  <table>\n    <thead><tr><th>Iter</th><th>x Approx</th><th>f(x)</th><th>Absolute Error</th><th>Relative Error</th><th>Notes</th></tr></thead>\n";
    out << "    <tbody>\n";
    for (const auto& rec : result.history) {
        out << "      <tr><td>" << rec.iteration << "</td><td>" << std::setprecision(8) << rec.x
            << "</td><td>" << rec.fx << "</td><td>" << std::scientific << rec.absoluteError
            << "</td><td>" << rec.relativeError << "</td><td>" << rec.note << "</td></tr>\n";
    }
    out << "    </tbody>\n  </table>\n";
    out << "</div>\n";

    out << "</body>\n</html>\n";
    return true;
}

} // namespace ees::ui

