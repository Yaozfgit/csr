#pragma once
#include <vector>
#include <cstddef>

class Matrix {
    public:
    //构造函数：n * n 方阵，初始为零
    explicit Matrix(std::size_t n) : n_(n), data_(n*n,0.0) {}

    //下标访问：读和写都可以用
    double& operator() (std::size_t i,std::size_t j);      //{return data_[i*n_+j];} 
    double operator() (std::size_t i,std::size_t j) const; //{return data_[i*n_+j];}   //const版本：只读

    //稀疏矩阵向量乘：y=A * x
    std::vector<double> multiply (const std::vector<double>& x) const;

    void print() const; //打印成易读的矩阵格式，用于调试
    
    double max() const;

    std::size_t size() const{ return n_; }   //内联小函数直接写在头文件里

    private:
    std::size_t n_;
    std::vector<double> data_;
};