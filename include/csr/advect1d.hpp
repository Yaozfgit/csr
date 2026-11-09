#pragma once
#include <vector>
#include <cstddef>
#include <csr/grid.hpp>

namespace csr {
    class LinearAdvection {
    public:
        LinearAdvection(const Grid1D& g, double a, const std::vector<double>& u0);
        void step(double dt);  //推进一个时间步
        const std::vector<double>& solution() const; // { return u_; }  //只读访问
        std::size_t n() const; // { return n_; }   //内部点数
    private:
        std::size_t n_;
        double a_;
        double dx_;
        std::vector<double> u_;   //内部边上的解，边界用周期条件补充
    };
}