#pragma once

#include <string>
#include <vector>
#include <map>

namespace ees::ui {

class CLIApp {
public:
    static int execute(int argc, char* argv[]);

private:
    static void printUsage();
    static void printVersion();

    // Command handlers
    static int handleRoot(const std::vector<std::string>& args);
    static int handleOptimize(const std::vector<std::string>& args);
    static int handleIntegrate(const std::vector<std::string>& args);
    static int handleInterpolate(const std::vector<std::string>& args);
    static int handleLinearSystem(const std::vector<std::string>& args);
    static int handleFactorize(const std::vector<std::string>& args);
    static int handleEigen(const std::vector<std::string>& args);
    static int handleSpline(const std::vector<std::string>& args);
    static int handleODE(const std::vector<std::string>& args);
    static int handlePDE(const std::vector<std::string>& args);
    static int handleCompare(const std::vector<std::string>& args);
    static int handleAdvise(const std::vector<std::string>& args);
    static int handleProblem(const std::vector<std::string>& args);
    static int handleTaylor(const std::vector<std::string>& args);

    // Parsing helpers
    static std::string getOption(const std::vector<std::string>& args, const std::string& opt, const std::string& defaultVal = "");
    static bool hasFlag(const std::vector<std::string>& args, const std::string& flag);
    static std::vector<double> parseCSV(const std::string& str);
    static std::vector<std::vector<double>> parseMatrix(const std::string& str);
};

} // namespace ees::ui

