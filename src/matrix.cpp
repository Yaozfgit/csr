#include "csr/matrix.hpp"

//成员变量初始化，成员n_初始化为n,data_初始化为n*n个0 
Matrix::Matrix(std::size_t n): n_(n), data_ (n * n, 0.0){}  

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

