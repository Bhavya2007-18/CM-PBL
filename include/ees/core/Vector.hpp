#pragma once

#include "Types.hpp"
#include <vector>
#include <initializer_list>
#include <iostream>
#include <string>
#include <stdexcept>

namespace ees::core {

class Vector {
public:
    Vector();
    explicit Vector(size_t size, double initVal = 0.0);
    Vector(std::initializer_list<double> list);
    Vector(const std::vector<double>& vec);
    Vector(std::vector<double>&& vec) noexcept;

    // Dimensions
    size_t size() const noexcept { return data_.size(); }
    bool empty() const noexcept { return data_.empty(); }
    void resize(size_t newSize, double val = 0.0) { data_.resize(newSize, val); }
    void clear() noexcept { data_.clear(); }

    // Raw pointers
    double* data() noexcept { return data_.data(); }
    const double* data() const noexcept { return data_.data(); }

    // Element access
    double& operator[](size_t index) noexcept { return data_[index]; }
    const double& operator[](size_t index) const noexcept { return data_[index]; }
    double& at(size_t index);
    const double& at(size_t index) const;

    // Iterators
    auto begin() noexcept { return data_.begin(); }
    auto end() noexcept { return data_.end(); }
    auto begin() const noexcept { return data_.begin(); }
    auto end() const noexcept { return data_.end(); }
    auto cbegin() const noexcept { return data_.cbegin(); }
    auto cend() const noexcept { return data_.cend(); }

    // Arithmetic operators
    Vector operator-() const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;
    Vector operator/(double scalar) const;

    Vector& operator+=(const Vector& other);
    Vector& operator-=(const Vector& other);
    Vector& operator*=(double scalar);
    Vector& operator/=(double scalar);

    bool operator==(const Vector& other) const;
    bool operator!=(const Vector& other) const;

    // Element-wise operations
    Vector cwiseProduct(const Vector& other) const;
    Vector cwiseQuotient(const Vector& other) const;

    // Vector operations
    double dot(const Vector& other) const;
    Vector cross(const Vector& other) const; // 3D cross product

    // Norms
    double norm1() const;       // Manhattan / L1
    double norm2() const;       // Euclidean / L2
    double norm() const { return norm2(); }
    double normInf() const;     // Infinity norm (max abs)
    Vector normalized() const;
    double distanceTo(const Vector& other) const;

    // Min / Max / Sum / Mean
    double minCoeff() const;
    double maxCoeff() const;
    double sum() const;
    double mean() const;

    // Generator helpers
    static Vector zero(size_t size);
    static Vector ones(size_t size);
    static Vector constant(size_t size, double value);
    static Vector unit(size_t size, size_t index);
    static Vector linspace(double start, double end, size_t numPoints);

    // Conversion & formatting
    const std::vector<double>& toStdVector() const noexcept { return data_; }
    std::string toString(int precision = 6) const;

    friend std::ostream& operator<<(std::ostream& os, const Vector& vec);
    friend Vector operator*(double scalar, const Vector& vec);

private:
    std::vector<double> data_;
};

Vector operator*(double scalar, const Vector& vec);

} // namespace ees::core
