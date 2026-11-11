#include <gtest/gtest.h>
#include <csr/grid.hpp>
#include <numeric>
#include <csr/advect1d.hpp>
#include <stdexcept>
#include <cmath>

double run_error(int n);
 // 网格加密一倍，误差应缩小约 4 倍（二阶）：误差比 ≈ (1/2)² = 1/4
double run_error(int n) {
    Grid1D g(n);
    std::cout << "n=" << n
              << " n_interior=" << g.n_interior()
              << " h=" << g.h()
              << " L=" << n * g.h() << "\n";
    std::vector<double> u0(g.n_interior());
    for (std::size_t i = 0; i < n; ++i)
        u0[i] = std::sin(2.0*M_PI*(i+0.5)*g.h()); //格子中心取值
    csr::LinearAdvection adv(g, 1.0, u0);
    const int nsteps = static_cast<int>(std::ceil(2.0*n)); //步数和n的关系
    //dt和h建立关系
    const double dt = 1.0/nsteps;
    for (int s = 0; s <nsteps; ++s) adv.step(dt);
    double err = 0;
    for (std::size_t i = 0; i<n; ++i)
        // TVD 格式 + L∞ = 结构性低估, L∞范式
        //err = std::max(err,std::abs(adv.solution()[i]-u0[i]));  
        err += std::abs(adv.solution()[i]-u0[i]); //L1范式
    return err * g.h();
}

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
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        num += (i + 0.5) * g.h() * u[i];   // grid_h = 1.0/200，或直接 g.h()
        den += u[i];
    }
    const double center = num / den;        // 波的质心，期望 0.5
    EXPECT_NEAR(center, 0.5, 3.0 / 200.0);
    /*
    std::size_t imax = 0;
    for (std::size_t i = 0; i<u.size(); ++i)
        if (u[i] > u[imax]) imax = i;
    EXPECT_NEAR(static_cast<double>(imax)/200.0, 0.5, 3.0 /200.0);

    //耗散：峰值低于1.0（物理，没有问题）
    EXPECT_LT(u[imax], 1.0);
    EXPECT_GT(u[imax], 0.5);   
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
 
TEST(AdvectionTest, SecondOrderConvergence) {   
    double e80 = run_error(80), e160 = run_error(160);
    double order = std::log2(e80/e160);
    std::cout<<"observed order = "<<order<<"\n";
    EXPECT_GE(order, 1.7); //理论是2.0，留点余量
}

TEST(AdvectionTest, NoOvershoot) {
    Grid1D g(100);
    std::vector<double> u0(g.n_interior(), 0.0);
    for (std::size_t i = 40; i < 60; ++i) u0[i]=1.0;
    csr::LinearAdvection adv(g, 1.0, u0);
    const int nsteps = 200;
    const double dt = 1.0/nsteps;
    for (int s = 0; s<nsteps; ++s) adv.step(dt);
    for (double v : adv.solution()) {
        EXPECT_LE(v, 1.0+1e-12);
        EXPECT_GE(v, -1e-12);
    }
}

TEST(AdvectionTest, WaveMovesRightward) {
    Grid1D g(100);
    std::vector<double> u0(g.n_interior(), 0.0);
    u0[30] = 1.0;                              // 高斯式单峰放在 x≈0.30
    csr::LinearAdvection adv(g, 1.0, u0);
    adv.step(0.005);                           // a=1，走 0.005 ≈ 半格
    const auto& u = adv.solution();
    // 波峰应出现在 30 或 31 号格（向 +x 移动），不许出现在 29
    EXPECT_GT(u[31] + u[30], u[29] + u[30] - 0.5*u[30]);  // 粗查：右侧比左侧保留了更多
}

TEST(GridTest, DomainLengthConsistency) {
    for (int n : {10, 80, 161}) {
        Grid1D g(n);
        EXPECT_DOUBLE_EQ(g.h() * g.n_interior(), 1.0);   // 域长恒为 1
    }
}
