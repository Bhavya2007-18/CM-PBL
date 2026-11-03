#pragma once

#include "Types.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>

namespace ees::core {

enum class TokenType {
    Number,
    Variable,
    Operator,
    Function,
    LeftParen,
    RightParen,
    Comma
};

struct Token {
    TokenType type;
    std::string text;
    double value{0.0};
};

class ExpressionParser {
public:
    ExpressionParser();
    explicit ExpressionParser(const std::string& expression);

    void setExpression(const std::string& expression);
    const std::string& getExpression() const noexcept { return rawExpr_; }

    // Variable management
    void setVariable(const std::string& name, double value);
    double getVariable(const std::string& name) const;

    // Evaluation
    double evaluate();
    double evaluate(double x);
    double evaluate(double x, double y);
    double evaluate(const std::unordered_map<std::string, double>& variables);
    double evaluate(const Vector& point, const std::vector<std::string>& varNames);

    // Differentiation utilities
    double derivativeCentral(double x, double h = 1e-6);
    double derivativeForward(double x, double h = 1e-6);
    double derivativeBackward(double x, double h = 1e-6);
    double secondDerivative(double x, double h = 1e-4);

    // Multivariate differentiation
    Vector gradient(const Vector& point, const std::vector<std::string>& varNames, double h = 1e-6);
    Matrix hessian(const Vector& point, const std::vector<std::string>& varNames, double h = 1e-4);

    // Safe evaluation check
    bool isValid() const noexcept { return isValid_; }
    const std::string& getErrorMessage() const noexcept { return errorMessage_; }

private:
    void tokenize();
    void parseToRPN();
    int getPrecedence(const std::string& op) const;
    bool isRightAssociative(const std::string& op) const;

    std::string rawExpr_;
    std::vector<Token> rpnQueue_;
    std::unordered_map<std::string, double> variables_;
    bool isValid_{false};
    std::string errorMessage_;
};

// Numerical differentiation helpers for arbitrary std::function
class NumericalDifferentiation {
public:
    static double centralDiff(const RealFunction& f, double x, double h = 1e-6);
    static double forwardDiff(const RealFunction& f, double x, double h = 1e-6);
    static double backwardDiff(const RealFunction& f, double x, double h = 1e-6);
    static double secondDiff(const RealFunction& f, double x, double h = 1e-4);
    static Vector gradient(const std::function<double(const Vector&)>& f, const Vector& x, double h = 1e-6);
    static Matrix hessian(const std::function<double(const Vector&)>& f, const Vector& x, double h = 1e-4);
};

} // namespace ees::core

