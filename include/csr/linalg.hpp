#pragma once
#include <csr/matrix.hpp>
#include <cstddef>
#include <vector>

//高斯消元（带部分主元），A是按值传入
std::vector<double> gauss_solve(Matrix A, std::vector<double> b); 