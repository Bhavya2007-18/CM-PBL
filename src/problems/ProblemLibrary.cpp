#include "ees/problems/ProblemLibrary.hpp"
#include "ees/unit1/RootFinding.hpp"
#include "ees/unit3/LinearSystems.hpp"
#include "ees/unit3/Splines.hpp"
#include "ees/unit4/PDESolvers.hpp"
#include <cmath>
#include <stdexcept>

namespace ees::problems {

std::vector<EngineeringProblem> ProblemLibrary::getAllProblems() {
    return {
        {
            "MECH-01",
            "Van der Waals Real Gas Molar Volume",
            EngineeringDomain::Mechanical,
            "Compute molar volume V of carbon dioxide at high pressure P and temperature T using the non-ideal Van der Waals cubic equation of state.",
            "(P + a / V^2) * (V - b) = R * T  ==>  P*V^3 - (P*b + R*T)*V^2 + a*V - a*b = 0",
            "P = 50 atm (5.066e6 Pa), T = 320 K, R = 8.314 J/(mol K), CO2 constants a = 0.3643 J m^3/mol^2, b = 4.267e-5 m^3/mol",
            "Newton-Raphson / Secant Method",
            "Molar Volume V ≈ 0.000378 m^3/mol"
        },
        {
            "CIVIL-01",
            "Pin-Jointed Planar Truss Force Equilibrium",
            EngineeringDomain::Civil,
            "Determine internal axial forces in structural members of a triangular truss supporting downward joint load.",
            "Sum(Fx) = 0, Sum(Fy) = 0 at free nodes  ==>  [A] * [F] = [Load]",
            "Load = 10 kN at joint, Member angles = 45 and 60 degrees",
            "Gauss Elimination with Partial Pivoting / LU Factorization",
            "Internal Forces F = [7.07 kN, -5.00 kN, 8.66 kN]"
        },
        {
            "ELEC-01",
            "Nonlinear Shockley Diode Circuit Operating Point",
            EngineeringDomain::Electrical,
            "Find the forward bias junction voltage V_D across a p-n diode in series with a resistor R and supply voltage V_s.",
            "f(V_D) = (V_s - V_D)/R - I_s * (exp(V_D / (n * V_t)) - 1) = 0",
            "V_s = 5.0 V, R = 1000 Ohm, I_s = 1e-12 A, n*V_t = 0.026 V",
            "Newton-Raphson (Quadratic Convergence)",
            "Operating Voltage V_D ≈ 0.684 V"
        },
        {
            "THERM-01",
            "Transient 1D Heat Conduction in Insulated Rod",
            EngineeringDomain::Thermal,
            "Model the time-dependent temperature profile in a cylindrical aluminum rod quenched in cold water at both ends.",
            "du/dt = alpha * d^2u/dx^2, u(0, t) = 0, u(L, t) = 0, u(x, 0) = 100 * sin(pi * x / L)",
            "Length L = 1.0 m, thermal diffusivity alpha = 9.7e-5 m^2/s, time t = 500 s",
            "Parabolic Finite Difference FTCS Scheme",
            "Centerline temperature diffuses smoothly from 100 C to ~62 C"
        },
        {
            "GEN-01",
            "Experimental Lubricating Oil Viscosity Calibration",
            EngineeringDomain::General,
            "Interpolate motor oil dynamic viscosity at intermediate temperature from experimental laboratory measurements.",
            "Natural Cubic Spline C^2 Piecewise Polynomial Approximation",
            "T = [0, 20, 40, 60, 80, 100] C, Viscosity mu = [1.79, 1.00, 0.65, 0.47, 0.36, 0.28] mPa.s",
            "Natural Cubic Spline / Gregory-Newton Backward",
            "Viscosity at T = 50 C ≈ 0.552 mPa.s"
        }
    };
}

std::vector<EngineeringProblem> ProblemLibrary::getProblemsByDomain(EngineeringDomain domain) {
    auto all = getAllProblems();
    std::vector<EngineeringProblem> filtered;
    for (const auto& p : all) {
        if (p.domain == domain) filtered.push_back(p);
    }
    return filtered;
}

EngineeringProblem ProblemLibrary::getProblemById(const std::string& id) {
    auto all = getAllProblems();
    for (const auto& p : all) {
        if (p.id == id) return p;
    }
    throw std::invalid_argument("Problem ID not found: " + id);
}

double ProblemLibrary::solveVanDerWaals(double P, double T, double R, double a, double b) {
    // Equation: f(V) = (P + a/(V^2))*(V - b) - R*T = 0
    auto f = [P, T, R, a, b](double V) {
        return (P + a / (V * V)) * (V - b) - R * T;
    };
    // Initial guess from ideal gas law V0 = R*T / P
    double vIdeal = (R * T) / P;
    auto res = unit1::RootFindingLab::newtonRaphson(f, vIdeal, 1e-8, 100);
    return res.solution;
}

core::Vector ProblemLibrary::solveTrussNodalForces() {
    // 3-bar truss equilibrium equations:
    // [ cos(45)  cos(60)   0      ] [ F1 ]   [ 0  ]
    // [ sin(45)  0         1.0    ] [ F2 ] = [ 10 ]
    // [ 0        sin(60)  -cos(45)] [ F3 ]   [ 0  ]
    double c45 = std::cos(3.141592653589793 / 4.0);
    double s45 = std::sin(3.141592653589793 / 4.0);
    double c60 = std::cos(3.141592653589793 / 3.0);
    double s60 = std::sin(3.141592653589793 / 3.0);

    core::Matrix A{
        { c45,  c60,  0.0},
        { s45,  0.0,  1.0},
        { 0.0,  s60, -c45}
    };
    core::Vector b{0.0, 10.0, 0.0};

    auto sol = unit3::LinearAlgebraLab::gaussElimination(A, b);
    return sol.solution;
}

double ProblemLibrary::solveDiodeOperatingPoint(double V_supply, double R_resistor, double I_s, double V_t) {
    // f(V_D) = (V_supply - V_D) / R_resistor - I_s * (exp(V_D / V_t) - 1) = 0
    auto f = [V_supply, R_resistor, I_s, V_t](double Vd) {
        return (V_supply - Vd) / R_resistor - I_s * (std::exp(Vd / V_t) - 1.0);
    };
    auto res = unit1::RootFindingLab::newtonRaphson(f, 0.6, 1e-8, 100);
    return res.solution;
}

core::Matrix ProblemLibrary::solveTransientRodCooling(double alpha, double L, double tFinal) {
    auto u0 = [](double x){ return 100.0 * std::sin(3.141592653589793 * x); };
    auto sol = unit4::PDESolverLab::solveHeatEquation1D(alpha, L, tFinal, 21, 51, u0, 0.0, 0.0);
    return sol.grid;
}

double ProblemLibrary::interpolateOilViscosity(double temperatureC) {
    core::Vector T{0.0, 20.0, 40.0, 60.0, 80.0, 100.0};
    core::Vector mu{1.79, 1.00, 0.65, 0.47, 0.36, 0.28};
    auto res = unit3::SplineLab::naturalCubicSpline(T, mu, temperatureC);
    return res.interpolatedY;
}

} // namespace ees::problems

