#pragma once

#include "ees/core/Types.hpp"
#include <string>

namespace ees::ui {

class Workbench {
public:
    Workbench();
    void run();

    // Lab Modules
    void showDashboard();
    void runUnit1Foundations();
    void runUnit1Roots();
    void runUnit1Optimization();

    void runUnit2Interpolation();
    void runUnit2Integration();

    void runUnit3LinearSystems();
    void runUnit3Eigenvalues();
    void runUnit3Splines();

    void runUnit4ODE();
    void runUnit4PDE();

    void runProblemLibrary();
    void runMethodComparison();
    void runMethodAdvisor();
    void runDemoMode(bool pauseAtEnd = true); // One-click comprehensive PBL Viva demonstration!

private:
    void printHeader(const std::string& title);
    void pause();
    int readInt(int minVal, int maxVal);
    double readDouble(const std::string& prompt, double defaultVal);
    std::string readString(const std::string& prompt, const std::string& defaultVal);
};

} // namespace ees::ui
