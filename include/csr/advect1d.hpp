#pragma once
#include <vector>
#include <cstddef>
#include <csr/grid.hpp>

namespace csr {
    class LinearAdvection {
    public:
        LinearAdvection(const Grid1D& g, double a, const std::vector<double>& u0);
        void step(double dt);  //推进一个时间步
        cosnt std::vector<double>& solution() const;  //只读访问
        std::size_t n() const;   //内部点数
    private:
        std::size_t n_;
        double a_;
        double dx_;
        std::vector<double> u_;   //内部垫上的解，边界用周期条件补充
    };
}