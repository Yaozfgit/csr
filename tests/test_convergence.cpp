#include <gtest/gtest.h>
#include <csr/poisson.hpp>
#include <csr/linalg.hpp>
#include <cmath>
#include <vector>

namespace {
    //解poisson方程并返回最大误差 max [u_h - sin(pi*x)]
    double poisson_max_error(std::size_t n) {
        CSRMatrix A_csr = build_poisson(n);
        Matrix A_dense(n);
        for (std::size_t i = 0; i<n; ++i)
            for (std::size_t j = 0; j<n; ++j)
                A_dense(i,j) = A_csr.at(i,j);

        std::vector<double> b(n);
        for (std::size_t i = 0; i<n; ++i) {
            double x = static_cast<double>(i+1)/(n+1);
            b[i] = M_PI*M_PI*sin(M_PI*x);
    }
    std::vector<double> u = gauss_solve(A_dense, b);
    
    double max_err = 0.0;
    for (std::size_t i = 0; i<n; ++i){
        double x = static_cast<double>(i+1)/(n+1);
        max_err = std::max(max_err, std::abs(u[i] - std::sin(M_PI*x)));
    }
    return max_err;
}
}

//C1:单点误差量级合理
TEST(ConvergenceTest, ErrorMagnitude) {
    double e100 = poisson_max_error(100);
    EXPECT_LT(e100, 1e-3);
}

//C2：核心测试——实测收敛阶数为2左右
TEST(ConvergenceTest,SecondOrderConvergence) {
    std::size_t n = 16;
    double prev = poisson_max_error(n);
    for (int level = 0; level<4; ++level) {
        n *=2;
        double curr = poisson_max_error(n);
        double p = std::log2(prev/curr);
        EXPECT_NEAR(p, 2.0,1.0);
        prev = curr;
    }
}

//C3:观察数据——打印表格，肉眼看收敛（顺手输出，不是断言）
TEST(ConvergenceTest,ConvergenceTable) {
    std::size_t n = 8;
    double prev = poisson_max_error(n);
    for (int level = 0; level < 4; ++level) {
        n *= 2;
        double curr =poisson_max_error(n);
        printf("n=%5zu err=%.3e order=%.3f\n", n, curr, std::log2(prev/curr));
        prev = curr;
    }
    SUCCEED();
}