#include <gtest/gtest.h>
#include <csr/grid.hpp>
#include <numeric>
#include <csr/advect1d.hpp>
#include <stdexcept>

double total(const std::vector<double>& u) {
    return std::accumulate(u.begin(), u.end(), 0.0);
}

// A1: 守恒性——周期边界下 Σu 严格不变（用相对量，容差机器精度级）
TEST(AdvectionTest, MassConserved) {
    Grid1D g(100);
    std::vector<double> u0(g.n_interior());
    for (std::size_t i = 0; i < u0.size(); ++i) u0[i] = 1.0 + 0.5 * i;
    csr::LinearAdvection adv(g, 1.0,u0);

    const double m0 = total(u0);
    for (int s = 0; s<200; ++s) adv.step(0.9*g.h() / 1.0);

    EXPECT_NEAR(total(adv.solution())/m0, 1.0, 1e-12);
}

// A2: 平移正确性——方波跑一个整周期后回到原位（周期域上 a=1、t=1 恰好绕一圈）
TEST(AdvectionTest, FullPeriodReturn) {
    // n=200，初始方波 [0.4,0.6]，t=1.0（恰好平移一个域长）
    Grid1D g(200);
    std::vector<double> u0(g.n_interior());
    //[80,120]即[0.4,0.6]之间u0=1
    for (std::size_t i = 80; i < 120; ++i) u0[i] = 1.0; 

    csr::LinearAdvection adv(g, 1.0, u0);
    const int nsteps = 250;
    const double dt = 1.0/nsteps;
    for (int s = 0; s < nsteps; ++s) adv.step(dt);

    const auto& u = adv.solution();

    //质量守恒（周期域，质量必须纹丝不动）
    // 精确解回到初始状态；Rusanov 有耗散，波形会变钝但**峰值位置不变**
    double m0 = 0.0, m1 = 0.0;
    for (double v : u0) m0 += v;
    for (double v : u) m1 += v;
    EXPECT_NEAR(m1/m0, 1.0, 1e-12);    
    // 断言：方波中心仍在 0.5 附近（±一个格），质量守恒，最大值略降（耗散的可视证据）
    std::size_t imax = 0;
    for (std::size_t i = 0; i<u.size(); ++i)
        if (u[i] > u[imax]) imax = i;
    EXPECT_NEAR(static_cast<double>(imax)/200.0, 0.5, 3.0 /200.0);

    //耗散：峰值低于1.0（物理，没有问题）
    EXPECT_LT(u[imax], 1.0);
    EXPECT_GT(u[imax], 0.5);

    /*/ 每 50 步打印一次方波轮廓（每 10 个点采样一个值）
    for (int s = 1; s <= 250; ++s) {
        adv.step(dt);
        if (s % 50 == 0) {
        std::cout << "t=" << dt * s << " : ";
        for (std::size_t i = 0; i < u.size(); i += 10) std::cout << u[i] << " ";
        std::cout << "\n";
    }
}
*/
}

// A3: CFL 违反 → 应该爆炸或被拒绝
TEST(AdvectionTest, CFLGuard) {
    // 决策题：step 里检测 dt > dx/|a| 时 throw，还是让用户自己负责？
    Grid1D g(100);
    std::vector<double> u0(g.n_interior(), 1.0);
    csr::LinearAdvection adv(g, 1.0, u0);

    EXPECT_THROW(adv.step(1.0), std::invalid_argument);   // dt=1.0 >> h=0.01 → 必炸
    EXPECT_NO_THROW(adv.step(0.009));                     // CFL=0.9 → 合法
    // 建议先写"手动大 dt 会得到 |u|→∞"的观察测试，再决定要不要加防护
}
    