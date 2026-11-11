#include <csr/advect1d.hpp>
#include <csr/flux.hpp>
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <csr/grid.hpp>
#include <iostream>

namespace csr {
LinearAdvection::LinearAdvection(const Grid1D& g, double a,
                                 const std::vector<double>& u0)
    : grid_(g), n_(g.n_interior()), a_(a), dx_(g.h()), u_(u0) {
        if (u0.size() != n_)
            throw std::invalid_argument("u0 size must match interior point");
    }
void csr::LinearAdvection::step(double dt) {
    if (dt > grid_.h() / std::abs(a_))
        throw std::invalid_argument("CFL violated: dt exceeds dx/|a|");

    const std::size_t n = u_.size();
    const auto rhs1 = compute_rhs(u_);              // 第一次：在 uⁿ 处算
    std::vector<double> u_star(n);
    for (std::size_t i = 0; i < n; ++i)
        u_star[i] = u_[i] - dt * rhs1[i];           // 预测：前进一整步

    const auto rhs2 = compute_rhs(u_star);          // 第二次：在 u* 处算！
    for (std::size_t i = 0; i < n; ++i)
        u_[i] = 0.5 * u_[i] + 0.5 * (u_star[i] - dt * rhs2[i]);  // 校正：平均
}

const std::vector<double>& LinearAdvection::solution() const { return u_; }
std::size_t LinearAdvection::n() const { return n_; }

static double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;                    // 异号（含零）→ 0
    return (std::abs(a) < std::abs(b)) ? a : b;      // 同号取小的
}
std::vector<double> csr::LinearAdvection::slopes(const std::vector<double>& u) const{          
    //const auto& u = u_;  //当前解
    const std::size_t n = u.size();
    const double h = grid_.h();
    std::vector<double> s(n,0.0);
    for (std::size_t i = 0; i <n; ++i) {
        std::size_t im = (i+n-1)%n;  //左边相邻，周期回绕
        std::size_t ip = (i+1)%n;   //右边相邻，周期回绕        
        s[i] = minmod(u[i]-u[im], u[ip]-u[i]) / h;
    }
    return s;
}
// 成员函数：std::vector<double> compute_rhs(const std::vector<double>& u) const;
std::vector<double> csr::LinearAdvection::compute_rhs(const std::vector<double>& u) const {
    const std::size_t n = u.size();
    const double h = grid_.h();
    // 1) 对"参数给的解"算斜率——注意是对 u 算，不是对 u_ 算！
    auto s = slopes(u);       
    // 2) 扫所有界面：算通量
    std::vector<double> flux(n);   // flux[i] 存界面 i 与 i+1 之间的 F_{i+½}
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t ip = (i + 1) % n;
        double uL = u[i]  + s[i]  * (0.5 * h);   // ← 改动②的外推公式，逐字对应
        double uR = u[ip] - s[ip] * (0.5 * h);
        flux[i] = rusanov_flux(uL, uR, a_);       // 通量函数零改动——F组测试仍在守护它
    }
    // 3) 散度：右界面通量 − 左界面通量
    std::vector<double> rhs(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t im = (i + n - 1) % n;
        rhs[i] = (flux[i] - flux[im]) / h;
    }
    return rhs;
}
}