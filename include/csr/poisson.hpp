#pragma once
#include <csr/csr_matrix.hpp> //头文件最小依赖是工程习惯，编译时间和大项目健康度的保证
#include <csr/grid.hpp>

CSRMatrix build_poisson(const Grid1D& g);