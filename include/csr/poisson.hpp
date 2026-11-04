#pragma once
#include <csr/csr_matrix.hpp> //头文件最小依赖是工程吸管，编译时间和大项目健康度的保证
#include <cstddef>

CSRMatrix build_poisson(std::size_t n);