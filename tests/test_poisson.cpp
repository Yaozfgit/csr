#include <gtest/gtest.h>
#include <csr/poisson.hpp>
#include <csr/linalg.hpp>
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

// B1: 小系统精确解（手算）
TEST(LinalgTest, SolveSmall) {
    Matrix A(2);
    A(0,0) = 2.0; A(0,1) = 1.0;
    A(1,0) = 1.0; A(1,1) = 3.0;
    auto x =gauss_solve(A, std::vector<double>{3.0,5.0});
    //手算:2x+y=3,x+3y=5;求解得：x=4/5,y=7/5
    EXPECT_DOUBLE_EQ(x[0],0.8);
    EXPECT_DOUBLE_EQ(x[1],1.4);
}

//B2：需要换主元的系统——对角元为0，不选主元被除数可能是0
TEST(LinalgTest, PivotingRequired) {
    Matrix A(2);
    A(0,0) = 0.0; A(0,1) = 1.0;  //第一行对角元是0 ！！
    A(1,0) = 1.0; A(1,1) = 1.0;
    auto x =gauss_solve(A, std::vector<double>{2.0,3.0});
    //手算:y=2, x=1
    EXPECT_DOUBLE_EQ(x[0],1.0);
    EXPECT_DOUBLE_EQ(x[1],2.0);
}

//B3: 验残差——解出来带回去，Ax-b==0(数值解的通用验收标准)
TEST(LinalgTest, ResidualSmall) {
    Matrix A(3);
    A(0,0) = 4; A(0,1) = 1; A(0,2) = 0;
    A(1,0) = 1; A(1,1) = 3; A(1,2) = 1;
    A(2,0) = 0; A(2,1) = 1; A(2,2) = 2;
    std::vector<double> b{1.0,2.0,3.0};
    auto x = gauss_solve(A, b);
    for (int i = 0; i < 3; ++i) {
        double r = 0;
        for (int j = 0; j<3; ++j) r +=A(i,j)*x[j];
        EXPECT_NEAR(r-b[i], 0.0, 1e-10);
    }
}