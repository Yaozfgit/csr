#include <gtest/gtest.h>
#include <csr/grid.hpp>
#include <numeric>
#include <csr/advect1d.hpp>

double total(const std::vector<double>& u) {
    return std::accumulate(u.begin(), u.end(), 0.0);
}

TEST(AdvectionTest, MassConserved) {
    Grid1D g(100);
    std::vector<double> u0(g.n_interior());
    for (std::size_t i = 0; i < u0.size(); ++i) u0[i] = 1.0 + 0.5 * i;
    csr::LinearAdvection adv(g, 1.0,u0);

    const double m0 = total(u0);
    for (int s = 0; s<200; ++s) adv.step(0.9*g.h() / 1.0);

    EXPECT_NEAR(total(adv.solution())/m0, 1.0, 1e-12);
}
/*
// A2: 平移正确性——方波跑一个整周期后回到原位（周期域上 a=1、t=1 恰好绕一圈）
TEST(AdvectionTest, FullPeriodReturn) {
    // n=200，初始方波 [0.4,0.6]，t=1.0（恰好平移一个域长）
    // 精确解回到初始状态；Rusanov 有耗散，波形会变钝但**峰值位置不变**
    // 断言：方波中心仍在 0.5 附近（±一个格），质量守恒，最大值略降（耗散的可视证据）
}

// A3: CFL 违反 → 应该爆炸或被拒绝
TEST(AdvectionTest, CFLGuard) {
    // 决策题：step 里检测 dt > dx/|a| 时 throw，还是让用户自己负责？
    // 建议先写"手动大 dt 会得到 |u|→∞"的观察测试，再决定要不要加防护
}
    */