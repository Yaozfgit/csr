#pragma once
#include <csr/grid.hpp>
#include <csr/csr_matrix.hpp>
#include <vector>

class BurgersSolver {
public:
    BurgersSolver(const Grid1D& g, double nu);

    // 扩散算子：Laplacian 矩阵（三对角，边界行带 ghost 外推）
    CSRMatrix build_diffusion() const;

    // 对流通量的离散右端项：-div(u^2/2)，显式
    // 界面通量用迎风平均（第一步先简单，正确优先）
    std::vector<double> convection_rhs(const std::vector<double>& u) const;

    // 推进一步：对流显式 + 扩散隐式
    // u_new = (I - dt*nu*L)^{-1} (u_old - dt*C(u_old))
    std::vector<double> step(const std::vector<double>& u, double dt) const;

    double h() const { return g_.h(); }

private:
    Grid1D g_;
    double nu_;
    std::size_t n_;
};