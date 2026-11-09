#include <csr/advect1d.hpp>
#include <csr/flux.hpp>
#include <stdexcept>
#include <cstddef>
#include <cmath>

namespace csr {
LinearAdvection::LinearAdvection(const Grid1D& g, double a,
                                 const std::vector<double>& u0)
    : n_(g.n_interior()), a_(a), dx_(g.h()), u_(u0) {
        if (u0.size() != n_)
            throw std::invalid_argument("u0 size must match interior point");
    }

void LinearAdvection::step(double dt) {
    if (dt > dx_ / std::abs(a_))
        throw std::invalid_argument("CFL violated: dt exceeds dx/|a|");
    // 第 1 阶段：算出所有 n 个界面通量（i 和 i+1 之间的界面记为 F[i]）
    std::vector<double> flux(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        std::size_t ip1 = (i + 1) % n_;          // 右邻居：周期回绕！
        flux[i] = rusanov_flux(u_[i], u_[ip1], a_);
    }
    // 第 2 阶段：统一更新（此时所有通量都是"旧解"算的）
    for (std::size_t i = 0; i < n_; ++i) {
        std::size_t im1 = (i + n_ - 1) % n_;     // 左邻居：i=0 时绕到 n-1
        u_[i] -= (dt / dx_) * (flux[i] - flux[im1]);
    }
}

const std::vector<double>& LinearAdvection::solution() const { return u_; }
std::size_t LinearAdvection::n() const { return n_; }
}