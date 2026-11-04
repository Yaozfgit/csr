#include <gtest/gtest.h>
#include <csr/poisson.hpp>
#include <csr/csr_matrix.hpp>
#include <cmath>

//---三对角结果————小n时逐元素核对组装结果--
TEST(PoissonTest, AssemblyTridiagonal) {
    const std::size_t n = 4;
    CSRMatrix A = build_poisson(n);
    //手算：h=1/5
    //对角线全为50，次对角线全部是-25，非零元数量nnz=2+3+3+2 =10
    EXPECT_EQ(A.nnz(), 10u);
    EXPECT_DOUBLE_EQ(A.at(0,0), 50.0);
    EXPECT_DOUBLE_EQ(A.at(0,1), -25.0);
    EXPECT_DOUBLE_EQ(A.at(2,1), -25.0);
    EXPECT_DOUBLE_EQ(A.at(3,3), 50.0);
    EXPECT_DOUBLE_EQ(A.at(0,3), 0.0);
}

//A2:差分算子作用于精确解解采样与f相同，解析解的对拍
TEST(PoissonTest, OperatorOnExactSolution) {
    const std::size_t n = 100;
    CSRMatrix A = build_poisson(n);
    std::vector<double> u(n), f(n);
    for (std::size_t i=0; i<n; ++i) {
        double x = static_cast<double>(i+1)/(n+1);
        u[i] = std::sin(M_PI * x);
        f[i] = M_PI * M_PI* std::sin(M_PI*x);
    }
    auto r = A.multiply(u);
    for (std::size_t i = 0; i<n; ++i)
        EXPECT_NEAR(r[i],f[i],0.1);
}