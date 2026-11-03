#include "ees/core/Vector.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

namespace ees::core {

Vector::Vector() = default;

Vector::Vector(size_t size, double initVal)
    : data_(size, initVal) {}

Vector::Vector(std::initializer_list<double> list)
    : data_(list) {}

Vector::Vector(const std::vector<double>& vec)
    : data_(vec) {}

Vector::Vector(std::vector<double>&& vec) noexcept
    : data_(std::move(vec)) {}

double& Vector::at(size_t index) {
    if (index >= data_.size()) {
        throw std::out_of_range("Vector index out of range: " + std::to_string(index) +
                                " >= " + std::to_string(data_.size()));
    }
    return data_[index];
}

const double& Vector::at(size_t index) const {
    if (index >= data_.size()) {
        throw std::out_of_range("Vector index out of range: " + std::to_string(index) +
                                " >= " + std::to_string(data_.size()));
    }
    return data_[index];
}

Vector Vector::operator-() const {
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        res[i] = -data_[i];
    }
    return res;
}

Vector Vector::operator+(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector dimensions do not match for addition: " +
                                    std::to_string(data_.size()) + " vs " + std::to_string(other.data_.size()));
    }
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        res[i] = data_[i] + other.data_[i];
    }
    return res;
}

Vector Vector::operator-(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector dimensions do not match for subtraction: " +
                                    std::to_string(data_.size()) + " vs " + std::to_string(other.data_.size()));
    }
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        res[i] = data_[i] - other.data_[i];
    }
    return res;
}

Vector Vector::operator*(double scalar) const {
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        res[i] = data_[i] * scalar;
    }
    return res;
}

Vector Vector::operator/(double scalar) const {
    if (std::abs(scalar) < 1e-30) {
        throw std::domain_error("Vector division by zero scalar");
    }
    double inv = 1.0 / scalar;
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        res[i] = data_[i] * inv;
    }
    return res;
}

Vector& Vector::operator+=(const Vector& other) {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector dimensions do not match for +=");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] += other.data_[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector& other) {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector dimensions do not match for -=");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] -= other.data_[i];
    }
    return *this;
}

Vector& Vector::operator*=(double scalar) {
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] *= scalar;
    }
    return *this;
}

Vector& Vector::operator/=(double scalar) {
    if (std::abs(scalar) < 1e-30) {
        throw std::domain_error("Vector division by zero scalar");
    }
    double inv = 1.0 / scalar;
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] *= inv;
    }
    return *this;
}

bool Vector::operator==(const Vector& other) const {
    if (data_.size() != other.data_.size()) return false;
    for (size_t i = 0; i < data_.size(); ++i) {
        if (std::abs(data_[i] - other.data_[i]) > 1e-12) return false;
    }
    return true;
}

bool Vector::operator!=(const Vector& other) const {
    return !(*this == other);
}

Vector Vector::cwiseProduct(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Dimension mismatch in cwiseProduct");
    }
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        res[i] = data_[i] * other.data_[i];
    }
    return res;
}

Vector Vector::cwiseQuotient(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Dimension mismatch in cwiseQuotient");
    }
    Vector res(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        if (std::abs(other.data_[i]) < 1e-30) {
            throw std::domain_error("Division by zero in cwiseQuotient");
        }
        res[i] = data_[i] / other.data_[i];
    }
    return res;
}

double Vector::dot(const Vector& other) const {
    if (data_.size() != other.data_.size()) {
        throw std::invalid_argument("Vector dimension mismatch in dot product: " +
                                    std::to_string(data_.size()) + " vs " + std::to_string(other.data_.size()));
    }
    double sum = 0.0;
    for (size_t i = 0; i < data_.size(); ++i) {
        sum += data_[i] * other.data_[i];
    }
    return sum;
}

Vector Vector::cross(const Vector& other) const {
    if (data_.size() != 3 || other.data_.size() != 3) {
        throw std::invalid_argument("Cross product is only defined for 3-dimensional vectors");
    }
    return Vector{
        data_[1] * other.data_[2] - data_[2] * other.data_[1],
        data_[2] * other.data_[0] - data_[0] * other.data_[2],
        data_[0] * other.data_[1] - data_[1] * other.data_[0]
    };
}

double Vector::norm1() const {
    double sum = 0.0;
    for (double val : data_) {
        sum += std::abs(val);
    }
    return sum;
}

double Vector::norm2() const {
    double sumSq = 0.0;
    for (double val : data_) {
        sumSq += val * val;
    }
    return std::sqrt(sumSq);
}

double Vector::normInf() const {
    double maxVal = 0.0;
    for (double val : data_) {
        maxVal = std::max(maxVal, std::abs(val));
    }
    return maxVal;
}

Vector Vector::normalized() const {
    double n = norm2();
    if (n < 1e-30) {
        throw std::runtime_error("Cannot normalize a zero vector");
    }
    return *this / n;
}

double Vector::distanceTo(const Vector& other) const {
    return (*this - other).norm2();
}

double Vector::minCoeff() const {
    if (data_.empty()) throw std::runtime_error("minCoeff called on empty vector");
    return *std::min_element(data_.begin(), data_.end());
}

double Vector::maxCoeff() const {
    if (data_.empty()) throw std::runtime_error("maxCoeff called on empty vector");
    return *std::max_element(data_.begin(), data_.end());
}

double Vector::sum() const {
    return std::accumulate(data_.begin(), data_.end(), 0.0);
}

double Vector::mean() const {
    if (data_.empty()) throw std::runtime_error("mean called on empty vector");
    return sum() / static_cast<double>(data_.size());
}

Vector Vector::zero(size_t size) {
    return Vector(size, 0.0);
}

Vector Vector::ones(size_t size) {
    return Vector(size, 1.0);
}

Vector Vector::constant(size_t size, double value) {
    return Vector(size, value);
}

Vector Vector::unit(size_t size, size_t index) {
    if (index >= size) {
        throw std::out_of_range("Unit vector index out of range");
    }
    Vector v(size, 0.0);
    v[index] = 1.0;
    return v;
}

Vector Vector::linspace(double start, double end, size_t numPoints) {
    if (numPoints == 0) return Vector();
    if (numPoints == 1) return Vector{ start };
    Vector res(numPoints);
    double step = (end - start) / static_cast<double>(numPoints - 1);
    for (size_t i = 0; i < numPoints; ++i) {
        res[i] = start + static_cast<double>(i) * step;
    }
    return res;
}

std::string Vector::toString(int precision) const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << "[";
    for (size_t i = 0; i < data_.size(); ++i) {
        ss << data_[i];
        if (i + 1 < data_.size()) ss << ", ";
    }
    ss << "]";
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const Vector& vec) {
    os << vec.toString();
    return os;
}

Vector operator*(double scalar, const Vector& vec) {
    return vec * scalar;
}

} // namespace ees::core
