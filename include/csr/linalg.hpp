#pragma once
#include <cstddef>
#include <vector>
#include <iosfwd>

class Matrix {
public:
    explicit Matrix(std::size_t n) : n_(n), data_(n*n,0.0) {}
    double& operator() (std::size_t i,std::size_t j)      {return data_[i*n_+j];} 
    double operator() (std::size_t i,std::size_t j) const {return data_[i*n_+j];}
    std::size_t size() const { return n_; }
private:
    std::size_t n_;
    std::vector<double> data_;
};

//高斯消元（带部分主元），A是按值传入
std::vector<double> gauss_solve(Matrix A, std::vector<double> b); 