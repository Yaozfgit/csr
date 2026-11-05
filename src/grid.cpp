#include <csr/grid.hpp>
#include <stdexcept>
#include <vector>

//---构造函数：数据源在这里确定---
Grid1D::Grid1D(std::size_t n_interior): n_(n_interior) {
    if (n_ == 0) throw std::invalid_argument("Grid1D: need at least 1 interior point");
}

//---以下是以n_为基础推导出来的---
std::size_t Grid1D::n_interior() const { return n_; }   //直接读n_
std::size_t Grid1D::n_total() const    { return n_+2; }   //基于n_的推导
double Grid1D::h() const                 { return 1.0/(n_+1); } //推导

double Grid1D::x(std::size_t i) const {
    if (i>n_+1) throw std::out_of_range("Grid1D::x: index beyond boundary");
    return i * h();   //x(i)=i*h；边界i=0,i=n+1自动为0和1的边界上
}

std::vector<double> Grid1D::interior_points() const {
    std::vector<double> xs;   //只读函数中建临时量数组，是合法的，不改变成员
    xs.reserve(n_);  //预留空间(预留n_个元素)，性能的小细节
    for (std::size_t i = 1; i <=n_; ++i) xs.push_back(x(i)); //在尾部添加元素，用于扩容
    return xs;
}