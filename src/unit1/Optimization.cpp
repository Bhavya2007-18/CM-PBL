#include "ees/unit1/Optimization.hpp"
#include "ees/core/ErrorAnalysis.hpp"
#include <cmath>
#include <chrono>
#include <algorithm>
#include <sstream>

namespace ees::unit1 {

core::NumericalResult OptimizationLab::fibonacciSearch(
    core::RealFunction f,
    double a,
    double b,
    double tol,
    int maxIter,
    core::OptimizationGoal goal
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Linear;
    result.convergenceRateExplanation = "Fibonacci search achieves optimal interval reduction ratio among all comparison-based 1D methods for a fixed number of function evaluations.";

    auto objective = [f, goal](double x) {
        return (goal == core::OptimizationGoal::Minimize) ? f(x) : -f(x);
    };

    // Generate Fibonacci sequence up to N such that (b - a) / F_N < tol
    std::vector<double> F{1.0, 1.0};
    while ((b - a) / F.back() > tol && static_cast<int>(F.size()) < maxIter) {
        F.push_back(F[F.size() - 1] + F[F.size() - 2]);
    }
    int N = static_cast<int>(F.size()) - 1;

    double L = b - a;
    double x1 = a + (F[N - 2] / F[N]) * L;
    double x2 = a + (F[N - 1] / F[N]) * L;
    double f1 = objective(x1);
    double f2 = objective(x2);

    for (int k = 1; k <= N - 2; ++k) {
        double intervalWidth = b - a;
        core::IterationRecord rec;
        rec.iteration = k;
        rec.x = (a + b) * 0.5;
        rec.fx = (goal == core::OptimizationGoal::Minimize) ? f(rec.x) : f(rec.x);
        rec.absoluteError = intervalWidth;
        rec.aux["a"] = a;
        rec.aux["b"] = b;
        rec.aux["x1"] = x1;
        rec.aux["x2"] = x2;
        rec.aux["interval_width"] = intervalWidth;
        result.addRecord(rec);

        if (f1 > f2) {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (F[N - k - 1] / F[N - k]) * (b - a);
            f2 = objective(x2);
        } else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (F[N - k - 2] / F[N - k]) * (b - a);
            f1 = objective(x1);
        }

        if (intervalWidth < tol) break;
    }

    result.success = true;
    result.solution = (a + b) * 0.5;
    result.status = core::ConvergenceStatus::Converged;
    result.finalError = b - a;
    result.message = "Fibonacci search converged after " + std::to_string(result.history.size()) + " evaluation cycles.";

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult OptimizationLab::goldenSectionSearch(
    core::RealFunction f,
    double a,
    double b,
    double tol,
    int maxIter,
    core::OptimizationGoal goal
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Linear;
    result.convergenceRateExplanation = "Golden Section Search guarantees interval shrinkage by factor (sqrt(5)-1)/2 ≈ 0.618034 at each step.";

    auto objective = [f, goal](double x) {
        return (goal == core::OptimizationGoal::Minimize) ? f(x) : -f(x);
    };

    const double phi = (std::sqrt(5.0) - 1.0) / 2.0; // 0.6180339887...
    double d = phi * (b - a);
    double x1 = b - d;
    double x2 = a + d;
    double f1 = objective(x1);
    double f2 = objective(x2);

    for (int iter = 1; iter <= maxIter; ++iter) {
        double intervalWidth = b - a;
        double xMid = 0.5 * (a + b);

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.x = xMid;
        rec.fx = f(xMid);
        rec.absoluteError = intervalWidth;
        rec.aux["a"] = a;
        rec.aux["b"] = b;
        rec.aux["x1"] = x1;
        rec.aux["x2"] = x2;
        rec.aux["width"] = intervalWidth;
        result.addRecord(rec);

        if (intervalWidth < tol) {
            result.success = true;
            result.solution = xMid;
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = intervalWidth;
            result.message = "Golden Section Search converged within tolerance.";
            break;
        }

        if (f1 > f2) {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = objective(x2);
        } else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = objective(x1);
        }
    }

    if (!result.success) {
        result.solution = 0.5 * (a + b);
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = b - a;
        result.message = "Maximum iterations reached before reaching target tolerance.";
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult OptimizationLab::newtonsMethod1D(
    core::RealFunction f,
    double x0,
    double tol,
    int maxIter
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Quadratic;
    result.convergenceRateExplanation = "Newton's 1D optimization finds roots of the gradient f'(x)=0 with quadratic local convergence.";

    double x = x0;

    for (int iter = 1; iter <= maxIter; ++iter) {
        double df = core::NumericalDifferentiation::centralDiff(f, x, 1e-5);
        double d2f = core::NumericalDifferentiation::secondDiff(f, x, 1e-4);

        if (std::abs(d2f) < 1e-15) {
            result.success = false;
            result.solution = x;
            result.status = core::ConvergenceStatus::DivisionByZero;
            result.message = "Second derivative vanishes (f''(x) ≈ 0) at iteration " + std::to_string(iter);
            break;
        }

        double step = df / d2f;
        double xNext = x - step;
        double err = std::abs(step);

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.x = xNext;
        rec.fx = f(xNext);
        rec.absoluteError = err;
        rec.aux["f'(x)"] = df;
        rec.aux["f''(x)"] = d2f;
        rec.aux["step"] = -step;
        result.addRecord(rec);

        if (err < tol || std::abs(df) < tol) {
            result.success = true;
            result.solution = xNext;
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = err;
            result.message = "Newton's optimization method converged after " + std::to_string(iter) + " iterations.";
            break;
        }

        x = xNext;
    }

    if (!result.success && result.status != core::ConvergenceStatus::DivisionByZero) {
        result.solution = x;
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = result.history.empty() ? 0.0 : result.history.back().absoluteError;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult OptimizationLab::steepestDescent(
    const std::function<double(const core::Vector&)>& f,
    const core::Vector& initialPoint,
    double stepSize,
    double tol,
    int maxIter,
    bool useBacktracking
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Linear;
    result.convergenceRateExplanation = "Steepest descent converges linearly, with rate governed by the condition number of the Hessian matrix.";

    core::Vector x = initialPoint;

    for (int iter = 1; iter <= maxIter; ++iter) {
        core::Vector grad = core::NumericalDifferentiation::gradient(f, x, 1e-5);
        double gradNorm = grad.norm2();

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.xVector = x;
        rec.x = x.size() > 0 ? x[0] : 0.0;
        rec.fx = f(x);
        rec.absoluteError = gradNorm;
        rec.aux["grad_norm"] = gradNorm;

        if (gradNorm < tol) {
            result.addRecord(rec);
            result.success = true;
            result.solutionVector = x;
            result.solution = rec.fx;
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = gradNorm;
            result.message = "Steepest Descent converged to stationary point (||grad|| < tol).";
            break;
        }

        // Search direction d = -grad
        core::Vector d = -grad;
        double alpha = stepSize;

        if (useBacktracking) {
            // Armijo condition backtracking line search: f(x + alpha*d) <= f(x) + c1 * alpha * grad^T * d
            double fCurrent = rec.fx;
            double c1 = 1e-4;
            double rho = 0.5;
            double directionalDeriv = grad.dot(d); // -||grad||^2

            for (int bt = 0; bt < 20; ++bt) {
                core::Vector xTrial = x + alpha * d;
                if (f(xTrial) <= fCurrent + c1 * alpha * directionalDeriv) {
                    break;
                }
                alpha *= rho;
            }
        }

        rec.aux["alpha"] = alpha;
        result.addRecord(rec);

        core::Vector xNext = x + alpha * d;
        if ((xNext - x).norm2() < 1e-14) {
            result.success = true;
            result.solutionVector = xNext;
            result.solution = f(xNext);
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = gradNorm;
            result.message = "Step size became smaller than machine precision.";
            break;
        }

        x = xNext;
    }

    if (!result.success) {
        result.solutionVector = x;
        result.solution = f(x);
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = result.history.empty() ? 0.0 : result.history.back().absoluteError;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

core::NumericalResult OptimizationLab::nelderMead(
    const std::function<double(const core::Vector&)>& f,
    const core::Vector& initialPoint,
    double initialStep,
    double tol,
    int maxIter
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    core::NumericalResult result;
    result.convergenceRate = core::ConvergenceRate::Sublinear;
    result.convergenceRateExplanation = "Nelder-Mead is a heuristic direct-search simplex method (gradient-free) with variable convergence.";

    const size_t n = initialPoint.size();
    if (n == 0) {
        result.success = false;
        result.status = core::ConvergenceStatus::InvalidInput;
        result.message = "Initial point has 0 dimensions.";
        return result;
    }

    // Initialize simplex with n + 1 vertices
    std::vector<core::Vector> simplex(n + 1, initialPoint);
    for (size_t i = 1; i <= n; ++i) {
        simplex[i][i - 1] += initialStep;
    }

    std::vector<double> fVals(n + 1);
    for (size_t i = 0; i <= n; ++i) {
        fVals[i] = f(simplex[i]);
    }

    // Nelder-Mead simplex transformation parameters
    const double alpha = 1.0; // Reflection
    const double gamma = 2.0; // Expansion
    const double rho   = 0.5; // Contraction
    const double sigma = 0.5; // Shrink

    for (int iter = 1; iter <= maxIter; ++iter) {
        // Sort simplex vertices so that f(x_0) <= f(x_1) <= ... <= f(x_n)
        std::vector<size_t> order(n + 1);
        for (size_t i = 0; i <= n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t i1, size_t i2) {
            return fVals[i1] < fVals[i2];
        });

        std::vector<core::Vector> sortedSimplex(n + 1);
        std::vector<double> sortedFVals(n + 1);
        for (size_t i = 0; i <= n; ++i) {
            sortedSimplex[i] = simplex[order[i]];
            sortedFVals[i] = fVals[order[i]];
        }
        simplex = sortedSimplex;
        fVals = sortedFVals;

        // Check both function value spread and simplex geometric size
        double fSpread = std::abs(fVals[n] - fVals[0]);
        double maxDist = 0.0;
        for (size_t i = 1; i <= n; ++i) {
            maxDist = std::max(maxDist, (simplex[i] - simplex[0]).norm2());
        }

        core::IterationRecord rec;
        rec.iteration = iter;
        rec.xVector = simplex[0];
        rec.x = simplex[0][0];
        rec.fx = fVals[0];
        rec.absoluteError = maxDist;
        rec.aux["best_f"] = fVals[0];
        rec.aux["worst_f"] = fVals[n];
        rec.aux["simplex_size"] = maxDist;
        rec.aux["f_spread"] = fSpread;

        if (maxDist < tol && fSpread < tol) {
            rec.note = "Converged (simplex collapsed & f spread < tol)";
            result.addRecord(rec);
            result.success = true;
            result.solutionVector = simplex[0];
            result.solution = fVals[0];
            result.status = core::ConvergenceStatus::Converged;
            result.finalError = maxDist;
            result.message = "Nelder-Mead simplex collapsed to optimal minimum within tolerance.";
            break;
        }

        // Calculate centroid of the best n points (excluding worst vertex x_n)
        core::Vector centroid = core::Vector::zero(n);
        for (size_t i = 0; i < n; ++i) {
            centroid += simplex[i];
        }
        centroid = centroid / static_cast<double>(n);

        // 1. Reflection: x_r = centroid + alpha * (centroid - x_n)
        core::Vector xr = centroid + alpha * (centroid - simplex[n]);
        double fxr = f(xr);

        if (fxr >= fVals[0] && fxr < fVals[n - 1]) {
            simplex[n] = xr;
            fVals[n] = fxr;
            rec.note = "Reflection";
        } else if (fxr < fVals[0]) {
            // 2. Expansion: x_e = centroid + gamma * (xr - centroid)
            core::Vector xe = centroid + gamma * (xr - centroid);
            double fxe = f(xe);
            if (fxe < fxr) {
                simplex[n] = xe;
                fVals[n] = fxe;
                rec.note = "Expansion";
            } else {
                simplex[n] = xr;
                fVals[n] = fxr;
                rec.note = "Reflection accepted";
            }
        } else {
            // 3. Contraction
            if (fxr < fVals[n]) {
                // Outside contraction
                core::Vector xc = centroid + rho * (xr - centroid);
                double fxc = f(xc);
                if (fxc <= fxr) {
                    simplex[n] = xc;
                    fVals[n] = fxc;
                    rec.note = "Outside Contraction";
                } else {
                    // Shrink
                    for (size_t i = 1; i <= n; ++i) {
                        simplex[i] = simplex[0] + sigma * (simplex[i] - simplex[0]);
                        fVals[i] = f(simplex[i]);
                    }
                    rec.note = "Shrink";
                }
            } else {
                // Inside contraction
                core::Vector xc = centroid - rho * (centroid - simplex[n]);
                double fxc = f(xc);
                if (fxc < fVals[n]) {
                    simplex[n] = xc;
                    fVals[n] = fxc;
                    rec.note = "Inside Contraction";
                } else {
                    // Shrink
                    for (size_t i = 1; i <= n; ++i) {
                        simplex[i] = simplex[0] + sigma * (simplex[i] - simplex[0]);
                        fVals[i] = f(simplex[i]);
                    }
                    rec.note = "Shrink";
                }
            }
        }

        result.addRecord(rec);
    }

    if (!result.success) {
        result.solutionVector = simplex[0];
        result.solution = fVals[0];
        result.status = core::ConvergenceStatus::MaxIterationsReached;
        result.finalError = result.history.empty() ? 0.0 : result.history.back().absoluteError;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    return result;
}

} // namespace ees::unit1
