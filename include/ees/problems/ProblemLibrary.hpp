#pragma once

#include "ees/core/Types.hpp"
#include "ees/core/Vector.hpp"
#include "ees/core/Matrix.hpp"
#include <string>
#include <vector>

namespace ees::problems {

enum class EngineeringDomain {
    Mechanical,
    Civil,
    Electrical,
    Thermal,
    General
};

struct EngineeringProblem {
    std::string id;
    std::string title;
    EngineeringDomain domain;
    std::string description;
    std::string mathematicalModel;
    std::string parametersExplanation;
    std::string recommendedMethods;
    std::string expectedResult;
};

class ProblemLibrary {
public:
    static std::vector<EngineeringProblem> getAllProblems();
    static std::vector<EngineeringProblem> getProblemsByDomain(EngineeringDomain domain);
    static EngineeringProblem getProblemById(const std::string& id);

    // Solvers for curated problems:
    // 1. Mechanical: Van der Waals molar volume solver
    static double solveVanDerWaals(double P, double T, double R, double a, double b);

    // 2. Civil: 3-bar truss nodal equilibrium system
    static core::Vector solveTrussNodalForces();

    // 3. Electrical: Shockley diode operating voltage in series circuit
    static double solveDiodeOperatingPoint(double V_supply, double R_resistor, double I_s, double V_t);

    // 4. Thermal: 1D Transient cooling rod temperature
    static core::Matrix solveTransientRodCooling(double alpha, double L, double tFinal);

    // 5. General: Experimental oil viscosity interpolation
    static double interpolateOilViscosity(double temperatureC);
};

} // namespace ees::problems

