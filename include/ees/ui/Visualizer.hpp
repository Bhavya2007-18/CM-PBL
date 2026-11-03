#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include "ees/core/NumericalResult.hpp"
#include <string>
#include <vector>

namespace ees::ui {

struct PlotSeries {
    std::string label;
    core::Vector x;
    core::Vector y;
    char pointChar{'*'};
    char lineChar{'-'};
};

class Visualizer {
public:
    // Render 2D function curve with axes, grid, and optional iteration marks
    static std::string plotFunction(
        core::RealFunction f,
        double xMin,
        double xMax,
        int width = 70,
        int height = 20,
        const std::vector<double>& highlightX = {},
        const std::string& title = "Function Curve"
    );

    // Plot Iteration vs Error (semi-log convergence graph)
    static std::string plotConvergence(
        const std::vector<core::IterationRecord>& history,
        int width = 70,
        int height = 18,
        const std::string& title = "Convergence Trajectory (Iteration vs Error)"
    );

    // Visualize shrinking search intervals (Bisection, Fibonacci, Golden Section)
    static std::string plotIntervalShrinkage(
        const std::vector<core::IterationRecord>& history,
        int barWidth = 60
    );

    // Multi-series comparative plot (e.g. comparing Euler vs RK4 vs Exact)
    static std::string plotMultiSeries(
        const std::vector<PlotSeries>& series,
        int width = 70,
        int height = 20,
        const std::string& title = "Trajectory Comparison"
    );

    // Export complete computation dashboard to standalone HTML5 with interactive Canvas/SVG charts
    static bool exportToHTML(
        const std::string& filePath,
        const std::string& pageTitle,
        const core::NumericalResult& result,
        const std::string& mathExplanation
    );
};

} // namespace ees::ui

