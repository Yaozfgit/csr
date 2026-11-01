#pragma once
#include <vector>
#include <cstddef>
#include <stdexcept>

class Matrix {
    public:
    //构造函数：n * n 方阵，初始为零
    explicit Matrix(std::size_t n);

    //下标访问：读和写都可以用
    double& operator()(std::size_t i, std::size_t j);
    double  operator()(std::size_t i, std::size_t j) const;   //const版本：只读

    //稀疏矩阵向量乘：y=A * x
    std::vector<double> multiply (const std::vector<double>& x) const;

    std::size_t size() const{ return n_; }   //内联小函数直接写在头文件里

    private:
    std::size_t n_;
    std::vector<double> data_;
};