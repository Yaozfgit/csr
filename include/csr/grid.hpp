#pragma once
#include <vector>
#include <cstddef>

//[0,1]上的均匀1D网格：内部点n个+两个边界点，合计n+2点
//左边界x(0)=0，x(i)=i*h (i=1..n 内部)，右边界 x(n+1)=1
class Grid1D {
public:
    explicit Grid1D(std::size_t n_interior);

    std::size_t n_interior() const;  //内部点数量n
    std::size_t n_total() const; //含边界点共n+2
    double h() const;  //网格间距1/(n+1)
    double x(std::size_t i) const; //第i个点的坐标（含边界：x(0)=0,x(n+1)=1)
    std::vector<double> interior_points() const; //内部点坐标数组，x_1...x_n

private:
    std::size_t n_;
};