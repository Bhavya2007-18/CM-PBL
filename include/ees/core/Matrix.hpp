#pragma once

#include "Types.hpp"
#include "Vector.hpp"
#include <vector>
#include <initializer_list>
#include <iostream>
#include <string>
#include <stdexcept>

namespace ees::core {

class Matrix {
public:
    Matrix();
    Matrix(size_t rows, size_t cols, double initVal = 0.0);
    Matrix(std::initializer_list<std::initializer_list<double>> list);
    Matrix(size_t rows, size_t cols, const std::vector<double>& data);

    // Dimensions
    size_t rows() const noexcept { return rows_; }
    size_t cols() const noexcept { return cols_; }
    bool isSquare() const noexcept { return rows_ == cols_ && rows_ > 0; }
    bool empty() const noexcept { return rows_ == 0 || cols_ == 0; }
    size_t size() const noexcept { return rows_ * cols_; }

    // Direct data access
    double* data() noexcept { return data_.data(); }
    const double* data() const noexcept { return data_.data(); }

    // Element access (0-indexed)
    double& operator()(size_t r, size_t c) noexcept { return data_[r * cols_ + c]; }
    const double& operator()(size_t r, size_t c) const noexcept { return data_[r * cols_ + c]; }
    double& at(size_t r, size_t c);
    const double& at(size_t r, size_t c) const;

    // Row and Column extraction/mutation
    Vector row(size_t r) const;
    Vector col(size_t c) const;
    void setRow(size_t r, const Vector& vec);
    void setCol(size_t c, const Vector& vec);

    // Arithmetic
    Matrix operator-() const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix operator/(double scalar) const;
    Matrix operator*(const Matrix& other) const;
    Vector operator*(const Vector& vec) const;

    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix& operator*=(double scalar);
    Matrix& operator/=(double scalar);

    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;

    // Matrix properties and transforms
    Matrix transpose() const;
    double trace() const;
    double determinant() const;

    // Elementary Row Operations (crucial for Gauss & Pivoting demonstration)
    void swapRows(size_t r1, size_t r2);
    void scaleRow(size_t r, double scalar);
    void addRowMultiple(size_t targetRow, size_t sourceRow, double factor); // target += factor * source

    // Submatrix & Augmentation
    Matrix augment(const Vector& colVec) const;
    Matrix augment(const Matrix& other) const;
    Matrix subMatrix(size_t startRow, size_t startCol, size_t numRows, size_t numCols) const;

    // Matrix Norms
    double normFrobenius() const;
    double normInf() const; // Max row sum
    double norm1() const;   // Max col sum

    // Factory Helpers
    static Matrix identity(size_t n);
    static Matrix zeros(size_t rows, size_t cols);
    static Matrix ones(size_t rows, size_t cols);
    static Matrix constant(size_t rows, size_t cols, double value);
    static Matrix diagonal(const Vector& diag);

    // String representation
    std::string toString(int precision = 6) const;
    std::string toAugmentedString(size_t separatorCol, int precision = 6) const;

    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);
    friend Matrix operator*(double scalar, const Matrix& mat);

private:
    size_t rows_{0};
    size_t cols_{0};
    std::vector<double> data_;
};

Matrix operator*(double scalar, const Matrix& mat);

} // namespace ees::core
