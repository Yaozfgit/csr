#include <csr/matrix.hpp>
#include <stdexcept>
#include <iostream>
#include <cstddef>

double& Matrix::operator()(std::size_t i, std::size_t j){
    if (i >= n_ || j>= n_) throw std::out_of_range("Matrix index out of range");
    return data_[i * n_ + j];
}

double Matrix::operator()(std::size_t i, std::size_t j) const {
    if (i >= n_ || j>= n_) throw std::out_of_range("Matrix index out of range");
    return data_[i * n_ + j];
}

std::vector<double> Matrix::multiply(const std::vector<double>& x) const {
    if (x.size() != n_) throw std::invalid_argument ("dimension mismatch");
    std::vector<double> y (n_, 0.0);
    for (std::size_t i = 0; i < n_; ++i)
    for (std::size_t j = 0; j < n_; ++j)
        y[i] += data_[i * n_ + j] *x[j];
    return y;
}

void Matrix::print() const {
    for (std::size_t i = 0; i<n_; ++i) {
        for (std::size_t j = 0; j<n_; ++j) {
            std::printf ("%8.3f ", data_[i*n_+j]);
        }
        std::printf("\n");
    }
}

double Matrix::max() const {
    double m = data_[0];
    for (std::size_t k = 1; k < data_.size(); ++k) {
        if (data_[k] >m) m=data_[k];
    }
    return m;
}