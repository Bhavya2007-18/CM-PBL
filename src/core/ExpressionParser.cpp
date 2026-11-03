#include "ees/core/ExpressionParser.hpp"
#include <cctype>
#include <cmath>
#include <stack>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace ees::core {

ExpressionParser::ExpressionParser() {
    variables_["pi"] = PI_VAL;
    variables_["e"] = E_VAL;
}

ExpressionParser::ExpressionParser(const std::string& expression)
    : ExpressionParser() {
    setExpression(expression);
}

void ExpressionParser::setExpression(const std::string& expression) {
    rawExpr_ = expression;
    try {
        tokenize();
        parseToRPN();
        isValid_ = true;
        errorMessage_.clear();
    } catch (const std::exception& ex) {
        isValid_ = false;
        errorMessage_ = ex.what();
    }
}

void ExpressionParser::setVariable(const std::string& name, double value) {
    variables_[name] = value;
}

double ExpressionParser::getVariable(const std::string& name) const {
    auto it = variables_.find(name);
    if (it != variables_.end()) {
        return it->second;
    }
    throw std::runtime_error("Undefined variable in ExpressionParser: " + name);
}

int ExpressionParser::getPrecedence(const std::string& op) const {
    if (op == "+" || op == "-") return 1;
    if (op == "*" || op == "/") return 2;
    if (op == "u-" || op == "u+") return 3; // Unary plus/minus
    if (op == "^") return 4;
    return 0;
}

bool ExpressionParser::isRightAssociative(const std::string& op) const {
    return op == "^" || op == "u-" || op == "u+";
}

void ExpressionParser::tokenize() {
    // We will produce tokens directly during parse or store intermediate token list
}

void ExpressionParser::parseToRPN() {
    rpnQueue_.clear();
    std::stack<Token> opStack;

    bool expectUnary = true;
    size_t i = 0;
    const size_t len = rawExpr_.length();

    while (i < len) {
        char c = rawExpr_[i];

        // Skip whitespace
        if (std::isspace(static_cast<unsigned char>(c))) {
            ++i;
            continue;
        }

        // Numbers (e.g. 12, 3.1415, 1e-5)
        if (std::isdigit(static_cast<unsigned char>(c)) || (c == '.' && i + 1 < len && std::isdigit(static_cast<unsigned char>(rawExpr_[i+1])))) {
            size_t start = i;
            bool dotSeen = false;
            bool expSeen = false;

            while (i < len) {
                char ch = rawExpr_[i];
                if (std::isdigit(static_cast<unsigned char>(ch))) {
                    ++i;
                } else if (ch == '.' && !dotSeen && !expSeen) {
                    dotSeen = true;
                    ++i;
                } else if ((ch == 'e' || ch == 'E') && !expSeen) {
                    expSeen = true;
                    ++i;
                    if (i < len && (rawExpr_[i] == '+' || rawExpr_[i] == '-')) {
                        ++i;
                    }
                } else {
                    break;
                }
            }
            std::string numStr = rawExpr_.substr(start, i - start);
            double val = std::stod(numStr);
            rpnQueue_.push_back(Token{TokenType::Number, numStr, val});
            expectUnary = false;
            continue;
        }

        // Identifiers: functions or variables or constants
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            size_t start = i;
            while (i < len && (std::isalnum(static_cast<unsigned char>(rawExpr_[i])) || rawExpr_[i] == '_')) {
                ++i;
            }
            std::string ident = rawExpr_.substr(start, i - start);
            std::string lowerIdent = ident;
            std::transform(lowerIdent.begin(), lowerIdent.end(), lowerIdent.begin(), [](unsigned char ch){ return std::tolower(ch); });

            // Check if function
            if (lowerIdent == "sin" || lowerIdent == "cos" || lowerIdent == "tan" ||
                lowerIdent == "asin" || lowerIdent == "acos" || lowerIdent == "atan" ||
                lowerIdent == "sinh" || lowerIdent == "cosh" || lowerIdent == "tanh" ||
                lowerIdent == "exp" || lowerIdent == "log" || lowerIdent == "ln" ||
                lowerIdent == "log10" || lowerIdent == "sqrt" || lowerIdent == "cbrt" ||
                lowerIdent == "abs" || lowerIdent == "floor" || lowerIdent == "ceil") {
                opStack.push(Token{TokenType::Function, lowerIdent, 0.0});
            } else if (lowerIdent == "pi") {
                rpnQueue_.push_back(Token{TokenType::Number, "pi", PI_VAL});
                expectUnary = false;
            } else if (lowerIdent == "e") {
                rpnQueue_.push_back(Token{TokenType::Number, "e", E_VAL});
                expectUnary = false;
            } else {
                rpnQueue_.push_back(Token{TokenType::Variable, ident, 0.0});
                expectUnary = false;
            }
            continue;
        }

        // Comma (for multi-arg functions)
        if (c == ',') {
            while (!opStack.empty() && opStack.top().type != TokenType::LeftParen) {
                rpnQueue_.push_back(opStack.top());
                opStack.pop();
            }
            if (opStack.empty()) {
                throw std::runtime_error("Mismatched parentheses or misplaced comma in expression");
            }
            expectUnary = true;
            ++i;
            continue;
        }

        // Left parenthesis
        if (c == '(') {
            opStack.push(Token{TokenType::LeftParen, "(", 0.0});
            expectUnary = true;
            ++i;
            continue;
        }

        // Right parenthesis
        if (c == ')') {
            while (!opStack.empty() && opStack.top().type != TokenType::LeftParen) {
                rpnQueue_.push_back(opStack.top());
                opStack.pop();
            }
            if (opStack.empty()) {
                throw std::runtime_error("Mismatched closing parenthesis ')' in expression");
            }
            opStack.pop(); // Pop '('

            // If the token at the top of the stack is a function token, pop it to RPN queue
            if (!opStack.empty() && opStack.top().type == TokenType::Function) {
                rpnQueue_.push_back(opStack.top());
                opStack.pop();
            }
            expectUnary = false;
            ++i;
            continue;
        }

        // Operators (+, -, *, /, ^)
        if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
            std::string op(1, c);
            if (expectUnary) {
                if (c == '-') op = "u-";
                else if (c == '+') op = "u+";
                else throw std::runtime_error("Unexpected operator '" + op + "'");
            }

            int prec = getPrecedence(op);
            bool rightAssoc = isRightAssociative(op);

            while (!opStack.empty() && opStack.top().type == TokenType::Operator) {
                int topPrec = getPrecedence(opStack.top().text);
                if ((!rightAssoc && prec <= topPrec) || (rightAssoc && prec < topPrec)) {
                    rpnQueue_.push_back(opStack.top());
                    opStack.pop();
                } else {
                    break;
                }
            }

            opStack.push(Token{TokenType::Operator, op, 0.0});
            expectUnary = true;
            ++i;
            continue;
        }

        throw std::runtime_error("Unexpected character in expression: '" + std::string(1, c) + "'");
    }

    while (!opStack.empty()) {
        if (opStack.top().type == TokenType::LeftParen) {
            throw std::runtime_error("Mismatched opening parenthesis '(' in expression");
        }
        rpnQueue_.push_back(opStack.top());
        opStack.pop();
    }
}

double ExpressionParser::evaluate() {
    if (!isValid_) {
        throw std::runtime_error("Cannot evaluate invalid expression: " + errorMessage_);
    }

    std::stack<double> valStack;

    for (const auto& token : rpnQueue_) {
        if (token.type == TokenType::Number) {
            valStack.push(token.value);
        } else if (token.type == TokenType::Variable) {
            valStack.push(getVariable(token.text));
        } else if (token.type == TokenType::Operator) {
            if (token.text == "u-") {
                if (valStack.empty()) throw std::runtime_error("Stack underflow on unary -");
                double a = valStack.top(); valStack.pop();
                valStack.push(-a);
            } else if (token.text == "u+") {
                // No change
            } else {
                if (valStack.size() < 2) throw std::runtime_error("Stack underflow on binary operator " + token.text);
                double b = valStack.top(); valStack.pop();
                double a = valStack.top(); valStack.pop();

                if (token.text == "+") {
                    valStack.push(a + b);
                } else if (token.text == "-") {
                    valStack.push(a - b);
                } else if (token.text == "*") {
                    valStack.push(a * b);
                } else if (token.text == "/") {
                    if (std::abs(b) < 1e-30) {
                        throw std::domain_error("Expression division by zero");
                    }
                    valStack.push(a / b);
                } else if (token.text == "^") {
                    if (a < 0.0 && std::floor(b) != b) {
                        throw std::domain_error("Negative base raised to non-integer power");
                    }
                    valStack.push(std::pow(a, b));
                }
            }
        } else if (token.type == TokenType::Function) {
            if (valStack.empty()) throw std::runtime_error("Stack underflow on function " + token.text);
            double a = valStack.top(); valStack.pop();

            if (token.text == "sin")   valStack.push(std::sin(a));
            else if (token.text == "cos")   valStack.push(std::cos(a));
            else if (token.text == "tan")   valStack.push(std::tan(a));
            else if (token.text == "asin") {
                if (a < -1.0 || a > 1.0) throw std::domain_error("asin out of domain [-1, 1]");
                valStack.push(std::asin(a));
            }
            else if (token.text == "acos") {
                if (a < -1.0 || a > 1.0) throw std::domain_error("acos out of domain [-1, 1]");
                valStack.push(std::acos(a));
            }
            else if (token.text == "atan")  valStack.push(std::atan(a));
            else if (token.text == "sinh")  valStack.push(std::sinh(a));
            else if (token.text == "cosh")  valStack.push(std::cosh(a));
            else if (token.text == "tanh")  valStack.push(std::tanh(a));
            else if (token.text == "exp")   valStack.push(std::exp(a));
            else if (token.text == "log" || token.text == "ln") {
                if (a <= 0.0) throw std::domain_error("logarithm of non-positive number: " + std::to_string(a));
                valStack.push(std::log(a));
            }
            else if (token.text == "log10") {
                if (a <= 0.0) throw std::domain_error("log10 of non-positive number: " + std::to_string(a));
                valStack.push(std::log10(a));
            }
            else if (token.text == "sqrt") {
                if (a < 0.0) throw std::domain_error("square root of negative number: " + std::to_string(a));
                valStack.push(std::sqrt(a));
            }
            else if (token.text == "cbrt")  valStack.push(std::cbrt(a));
            else if (token.text == "abs")   valStack.push(std::abs(a));
            else if (token.text == "floor") valStack.push(std::floor(a));
            else if (token.text == "ceil")  valStack.push(std::ceil(a));
            else throw std::runtime_error("Unknown function: " + token.text);
        }
    }

    if (valStack.size() != 1) {
        throw std::runtime_error("Invalid expression: multiple values remained on stack");
    }

    return valStack.top();
}

double ExpressionParser::evaluate(double x) {
    setVariable("x", x);
    return evaluate();
}

double ExpressionParser::evaluate(double x, double y) {
    setVariable("x", x);
    setVariable("y", y);
    return evaluate();
}

double ExpressionParser::evaluate(const std::unordered_map<std::string, double>& variables) {
    for (const auto& [k, v] : variables) {
        setVariable(k, v);
    }
    return evaluate();
}

double ExpressionParser::evaluate(const Vector& point, const std::vector<std::string>& varNames) {
    if (point.size() != varNames.size()) {
        throw std::invalid_argument("Vector size does not match variable names count");
    }
    for (size_t i = 0; i < point.size(); ++i) {
        setVariable(varNames[i], point[i]);
    }
    return evaluate();
}

double ExpressionParser::derivativeCentral(double x, double h) {
    double f1 = evaluate(x + h);
    double f0 = evaluate(x - h);
    return (f1 - f0) / (2.0 * h);
}

double ExpressionParser::derivativeForward(double x, double h) {
    double f1 = evaluate(x + h);
    double f0 = evaluate(x);
    return (f1 - f0) / h;
}

double ExpressionParser::derivativeBackward(double x, double h) {
    double f0 = evaluate(x);
    double f_1 = evaluate(x - h);
    return (f0 - f_1) / h;
}

double ExpressionParser::secondDerivative(double x, double h) {
    double f1 = evaluate(x + h);
    double f0 = evaluate(x);
    double f_1 = evaluate(x - h);
    return (f1 - 2.0 * f0 + f_1) / (h * h);
}

Vector ExpressionParser::gradient(const Vector& point, const std::vector<std::string>& varNames, double h) {
    const size_t n = point.size();
    Vector grad(n);
    Vector pForward = point;
    Vector pBackward = point;

    for (size_t i = 0; i < n; ++i) {
        pForward[i] = point[i] + h;
        pBackward[i] = point[i] - h;

        double fPlus = evaluate(pForward, varNames);
        double fMinus = evaluate(pBackward, varNames);

        grad[i] = (fPlus - fMinus) / (2.0 * h);

        pForward[i] = point[i];
        pBackward[i] = point[i];
    }
    return grad;
}

Matrix ExpressionParser::hessian(const Vector& point, const std::vector<std::string>& varNames, double h) {
    const size_t n = point.size();
    Matrix H(n, n);
    Vector p = point;

    double f0 = evaluate(point, varNames);

    for (size_t i = 0; i < n; ++i) {
        // Diagonal: d^2f / dx_i^2
        p[i] = point[i] + h;
        double fPlus = evaluate(p, varNames);
        p[i] = point[i] - h;
        double fMinus = evaluate(p, varNames);
        p[i] = point[i];

        H(i, i) = (fPlus - 2.0 * f0 + fMinus) / (h * h);

        // Off-diagonal: d^2f / (dx_i dx_j)
        for (size_t j = i + 1; j < n; ++j) {
            p[i] = point[i] + h; p[j] = point[j] + h;
            double f_pp = evaluate(p, varNames);

            p[i] = point[i] + h; p[j] = point[j] - h;
            double f_pm = evaluate(p, varNames);

            p[i] = point[i] - h; p[j] = point[j] + h;
            double f_mp = evaluate(p, varNames);

            p[i] = point[i] - h; p[j] = point[j] - h;
            double f_mm = evaluate(p, varNames);

            p[i] = point[i]; p[j] = point[j];

            double val = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h * h);
            H(i, j) = val;
            H(j, i) = val; // Symmetry by Clairaut's theorem
        }
    }
    return H;
}

// NumericalDifferentiation static implementations
double NumericalDifferentiation::centralDiff(const RealFunction& f, double x, double h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double NumericalDifferentiation::forwardDiff(const RealFunction& f, double x, double h) {
    return (f(x + h) - f(x)) / h;
}

double NumericalDifferentiation::backwardDiff(const RealFunction& f, double x, double h) {
    return (f(x) - f(x - h)) / h;
}

double NumericalDifferentiation::secondDiff(const RealFunction& f, double x, double h) {
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
}

Vector NumericalDifferentiation::gradient(const std::function<double(const Vector&)>& f, const Vector& x, double h) {
    const size_t n = x.size();
    Vector g(n);
    Vector pForward = x;
    Vector pBackward = x;

    for (size_t i = 0; i < n; ++i) {
        pForward[i] = x[i] + h;
        pBackward[i] = x[i] - h;
        g[i] = (f(pForward) - f(pBackward)) / (2.0 * h);
        pForward[i] = x[i];
        pBackward[i] = x[i];
    }
    return g;
}

Matrix NumericalDifferentiation::hessian(const std::function<double(const Vector&)>& f, const Vector& x, double h) {
    const size_t n = x.size();
    Matrix H(n, n);
    Vector p = x;
    double f0 = f(x);

    for (size_t i = 0; i < n; ++i) {
        p[i] = x[i] + h;
        double fPlus = f(p);
        p[i] = x[i] - h;
        double fMinus = f(p);
        p[i] = x[i];

        H(i, i) = (fPlus - 2.0 * f0 + fMinus) / (h * h);

        for (size_t j = i + 1; j < n; ++j) {
            p[i] = x[i] + h; p[j] = x[j] + h;
            double f_pp = f(p);

            p[i] = x[i] + h; p[j] = x[j] - h;
            double f_pm = f(p);

            p[i] = x[i] - h; p[j] = x[j] + h;
            double f_mp = f(p);

            p[i] = x[i] - h; p[j] = x[j] - h;
            double f_mm = f(p);

            p[i] = x[i]; p[j] = x[j];

            double val = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h * h);
            H(i, j) = val;
            H(j, i) = val;
        }
    }
    return H;
}

} // namespace ees::core

