#include "ees/core/Matrix.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace ees::core {

Matrix::Matrix() = default;

Matrix::Matrix(size_t rows, size_t cols, double initVal)
    : rows_(rows), cols_(cols), data_(rows * cols, initVal) {}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
    rows_ = list.size();
    cols_ = rows_ > 0 ? list.begin()->size() : 0;
    data_.reserve(rows_ * cols_);
    for (const auto& rowList : list) {
        if (rowList.size() != cols_) {
            throw std::invalid_argument("Inconsistent row sizes in matrix initializer list");
        }
        for (double val : rowList) {
            data_.push_back(val);
        }
    }
}

Matrix::Matrix(size_t rows, size_t cols, const std::vector<double>& data)
    : rows_(rows), cols_(cols), data_(data) {
    if (data_.size() != rows_ * cols_) {
        throw std::invalid_argument("Data size does not match rows * cols in Matrix construction");
    }
}

double& Matrix::at(size_t r, size_t c) {
    if (r >= rows_ || c >= cols_) {
        throw std::out_of_range("Matrix subscript out of range: (" +
                                std::to_string(r) + ", " + std::to_string(c) + ") in " +
                                std::to_string(rows_) + "x" + std::to_string(cols_));
    }
    return data_[r * cols_ + c];
}

const double& Matrix::at(size_t r, size_t c) const {
    if (r >= rows_ || c >= cols_) {
        throw std::out_of_range("Matrix subscript out of range: (" +
                                std::to_string(r) + ", " + std::to_string(c) + ") in " +
                                std::to_string(rows_) + "x" + std::to_string(cols_));
    }
    return data_[r * cols_ + c];
}

Vector Matrix::row(size_t r) const {
    if (r >= rows_) throw std::out_of_range("Row index out of range: " + std::to_string(r));
    Vector v(cols_);
    for (size_t c = 0; c < cols_; ++c) {
        v[c] = (*this)(r, c);
    }
    return v;
}

Vector Matrix::col(size_t c) const {
    if (c >= cols_) throw std::out_of_range("Col index out of range: " + std::to_string(c));
    Vector v(rows_);
    for (size_t r = 0; r < rows_; ++r) {
        v[r] = (*this)(r, c);
    }
    return v;
}

void Matrix::setRow(size_t r, const Vector& vec) {
    if (r >= rows_) throw std::out_of_range("Row index out of range: " + std::to_string(r));
    if (vec.size() != cols_) throw std::invalid_argument("Vector size does not match matrix columns");
    for (size_t c = 0; c < cols_; ++c) {
        (*this)(r, c) = vec[c];
    }
}

void Matrix::setCol(size_t c, const Vector& vec) {
    if (c >= cols_) throw std::out_of_range("Col index out of range: " + std::to_string(c));
    if (vec.size() != rows_) throw std::invalid_argument("Vector size does not match matrix rows");
    for (size_t r = 0; r < rows_; ++r) {
        (*this)(r, c) = vec[r];
    }
}

Matrix Matrix::operator-() const {
    Matrix res(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) {
        res.data_[i] = -data_[i];
    }
    return res;
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix dimension mismatch for addition");
    }
    Matrix res(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) {
        res.data_[i] = data_[i] + other.data_[i];
    }
    return res;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix dimension mismatch for subtraction");
    }
    Matrix res(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) {
        res.data_[i] = data_[i] - other.data_[i];
    }
    return res;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix res(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) {
        res.data_[i] = data_[i] * scalar;
    }
    return res;
}

Matrix Matrix::operator/(double scalar) const {
    if (std::abs(scalar) < 1e-30) throw std::domain_error("Matrix division by zero scalar");
    double inv = 1.0 / scalar;
    Matrix res(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) {
        res.data_[i] = data_[i] * inv;
    }
    return res;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("Matrix multiplication inner dimension mismatch: " +
                                    std::to_string(cols_) + " != " + std::to_string(other.rows_));
    }
    Matrix res(rows_, other.cols_, 0.0);
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t k = 0; k < cols_; ++k) {
            double rVal = (*this)(i, k);
            for (size_t j = 0; j < other.cols_; ++j) {
                res(i, j) += rVal * other(k, j);
            }
        }
    }
    return res;
}

Vector Matrix::operator*(const Vector& vec) const {
    if (cols_ != vec.size()) {
        throw std::invalid_argument("Matrix-Vector multiplication dimension mismatch: " +
                                    std::to_string(cols_) + " vs " + std::to_string(vec.size()));
    }
    Vector res(rows_, 0.0);
    for (size_t i = 0; i < rows_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < cols_; ++j) {
            sum += (*this)(i, j) * vec[j];
        }
        res[i] = sum;
    }
    return res;
}

Matrix& Matrix::operator+=(const Matrix& other) {
    *this = *this + other;
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& other) {
    *this = *this - other;
    return *this;
}

Matrix& Matrix::operator*=(double scalar) {
    *this = *this * scalar;
    return *this;
}

Matrix& Matrix::operator/=(double scalar) {
    *this = *this / scalar;
    return *this;
}

bool Matrix::operator==(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) return false;
    for (size_t i = 0; i < data_.size(); ++i) {
        if (std::abs(data_[i] - other.data_[i]) > 1e-12) return false;
    }
    return true;
}

bool Matrix::operator!=(const Matrix& other) const {
    return !(*this == other);
}

Matrix Matrix::transpose() const {
    Matrix res(cols_, rows_);
    for (size_t r = 0; r < rows_; ++r) {
        for (size_t c = 0; c < cols_; ++c) {
            res(c, r) = (*this)(r, c);
        }
    }
    return res;
}

double Matrix::trace() const {
    if (!isSquare()) throw std::domain_error("Trace is only defined for square matrices");
    double tr = 0.0;
    for (size_t i = 0; i < rows_; ++i) {
        tr += (*this)(i, i);
    }
    return tr;
}

double Matrix::determinant() const {
    if (!isSquare()) throw std::domain_error("Determinant is only defined for square matrices");
    const size_t n = rows_;
    if (n == 1) return (*this)(0, 0);
    if (n == 2) return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);

    // Compute determinant via LU / Gaussian Elimination with partial pivoting
    Matrix temp = *this;
    double det = 1.0;
    int sign = 1;

    for (size_t i = 0; i < n; ++i) {
        // Find pivot
        size_t maxRow = i;
        double maxVal = std::abs(temp(i, i));
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(temp(k, i)) > maxVal) {
                maxVal = std::abs(temp(k, i));
                maxRow = k;
            }
        }

        if (maxVal < 1e-14) {
            return 0.0; // Matrix is singular
        }

        if (maxRow != i) {
            temp.swapRows(i, maxRow);
            sign = -sign;
        }

        det *= temp(i, i);

        for (size_t k = i + 1; k < n; ++k) {
            double factor = temp(k, i) / temp(i, i);
            for (size_t j = i + 1; j < n; ++j) {
                temp(k, j) -= factor * temp(i, j);
            }
        }
    }

    return sign * det;
}

void Matrix::swapRows(size_t r1, size_t r2) {
    if (r1 >= rows_ || r2 >= rows_) throw std::out_of_range("Row index out of range for swap");
    if (r1 == r2) return;
    for (size_t c = 0; c < cols_; ++c) {
        std::swap((*this)(r1, c), (*this)(r2, c));
    }
}

void Matrix::scaleRow(size_t r, double scalar) {
    if (r >= rows_) throw std::out_of_range("Row index out of range for scaleRow");
    for (size_t c = 0; c < cols_; ++c) {
        (*this)(r, c) *= scalar;
    }
}

void Matrix::addRowMultiple(size_t targetRow, size_t sourceRow, double factor) {
    if (targetRow >= rows_ || sourceRow >= rows_) {
        throw std::out_of_range("Row index out of range in addRowMultiple");
    }
    for (size_t c = 0; c < cols_; ++c) {
        (*this)(targetRow, c) += factor * (*this)(sourceRow, c);
    }
}

Matrix Matrix::augment(const Vector& colVec) const {
    if (colVec.size() != rows_) {
        throw std::invalid_argument("Augmented vector length must match matrix row count");
    }
    Matrix res(rows_, cols_ + 1);
    for (size_t r = 0; r < rows_; ++r) {
        for (size_t c = 0; c < cols_; ++c) {
            res(r, c) = (*this)(r, c);
        }
        res(r, cols_) = colVec[r];
    }
    return res;
}

Matrix Matrix::augment(const Matrix& other) const {
    if (other.rows() != rows_) {
        throw std::invalid_argument("Augmented matrix row count must match");
    }
    Matrix res(rows_, cols_ + other.cols());
    for (size_t r = 0; r < rows_; ++r) {
        for (size_t c = 0; c < cols_; ++c) {
            res(r, c) = (*this)(r, c);
        }
        for (size_t c = 0; c < other.cols(); ++c) {
            res(r, cols_ + c) = other(r, c);
        }
    }
    return res;
}

Matrix Matrix::subMatrix(size_t startRow, size_t startCol, size_t numRows, size_t numCols) const {
    if (startRow + numRows > rows_ || startCol + numCols > cols_) {
        throw std::out_of_range("Submatrix boundary exceeds matrix dimensions");
    }
    Matrix res(numRows, numCols);
    for (size_t r = 0; r < numRows; ++r) {
        for (size_t c = 0; c < numCols; ++c) {
            res(r, c) = (*this)(startRow + r, startCol + c);
        }
    }
    return res;
}

double Matrix::normFrobenius() const {
    double sumSq = 0.0;
    for (double v : data_) {
        sumSq += v * v;
    }
    return std::sqrt(sumSq);
}

double Matrix::normInf() const {
    double maxRowSum = 0.0;
    for (size_t r = 0; r < rows_; ++r) {
        double rowSum = 0.0;
        for (size_t c = 0; c < cols_; ++c) {
            rowSum += std::abs((*this)(r, c));
        }
        maxRowSum = std::max(maxRowSum, rowSum);
    }
    return maxRowSum;
}

double Matrix::norm1() const {
    double maxColSum = 0.0;
    for (size_t c = 0; c < cols_; ++c) {
        double colSum = 0.0;
        for (size_t r = 0; r < rows_; ++r) {
            colSum += std::abs((*this)(r, c));
        }
        maxColSum = std::max(maxColSum, colSum);
    }
    return maxColSum;
}

Matrix Matrix::identity(size_t n) {
    Matrix res(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        res(i, i) = 1.0;
    }
    return res;
}

Matrix Matrix::zeros(size_t rows, size_t cols) {
    return Matrix(rows, cols, 0.0);
}

Matrix Matrix::ones(size_t rows, size_t cols) {
    return Matrix(rows, cols, 1.0);
}

Matrix Matrix::constant(size_t rows, size_t cols, double value) {
    return Matrix(rows, cols, value);
}

Matrix Matrix::diagonal(const Vector& diag) {
    size_t n = diag.size();
    Matrix res(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        res(i, i) = diag[i];
    }
    return res;
}

std::string Matrix::toString(int precision) const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision);
    for (size_t r = 0; r < rows_; ++r) {
        ss << "[ ";
        for (size_t c = 0; c < cols_; ++c) {
            ss << std::setw(precision + 4) << (*this)(r, c);
            if (c + 1 < cols_) ss << " ";
        }
        ss << " ]";
        if (r + 1 < rows_) ss << "\n";
    }
    return ss.str();
}

std::string Matrix::toAugmentedString(size_t separatorCol, int precision) const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision);
    for (size_t r = 0; r < rows_; ++r) {
        ss << "[ ";
        for (size_t c = 0; c < cols_; ++c) {
            if (c == separatorCol) {
                ss << "| ";
            }
            ss << std::setw(precision + 4) << (*this)(r, c);
            if (c + 1 < cols_) ss << " ";
        }
        ss << " ]";
        if (r + 1 < rows_) ss << "\n";
    }
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
    os << mat.toString();
    return os;
}

Matrix operator*(double scalar, const Matrix& mat) {
    return mat * scalar;
}

} // namespace ees::core
