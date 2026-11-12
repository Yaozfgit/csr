#include <gtest/gtest.h>
#include <csr/thomas.hpp>
#include <vector>
#include <cmath>
#include <iostream>

namespace {

std::vector<double> assemble_b(const std::vector<double>& a,
                               const std::vector<double>& d,
                               const std::vector<double>& c,
                               const std::vector<double>& x_exact) {
    const std::size_t n = d.size();
    std::vector<double> b(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        b[i] = d[i] * x_exact[i];
        if (i > 0)     b[i] += a[i] * x_exact[i-1];
        if (i + 1 < n) b[i] += c[i] * x_exact[i+1];
    }
    return b;
}

}  // namespace

TEST(ThomasTest, Tiny2x2) {
    const std::size_t n = 2;
    std::vector<double> a(n, -1.0), d(n, 4.0), c(n, -1.0);
    a[0] = 0.0; c[n-1] = 0.0;                    // ← 只清这两处
    std::vector<double> x_exact{1.0, 2.0};
    auto b = assemble_b(a, d, c, x_exact);       // 期望 b = {2, 7}
    EXPECT_NEAR(b[0], 2.0, 1e-14);               // 组装自检，错在这层就别怪求解器
    EXPECT_NEAR(b[1], 7.0, 1e-14);

    auto x = thomas_solve(a, d, c, b);
    EXPECT_NEAR(x[0], 1.0, 1e-12);
    EXPECT_NEAR(x[1], 2.0, 1e-12);
}

TEST(ThomasTest, MatchesDenseSolver) {
    const std::size_t n = 2;
    std::vector<double> a(n, -1.0), d(n, 4.0), c(n, -1.0);
    a[0] = 0.0; c[n-1] = 0.0;
    std::vector<double> x_exact(n);
    for (std::size_t i = 0; i < n; ++i)
        x_exact[i] = std::sin(0.7 * i) + 0.3;
    auto b = assemble_b(a, d, c, x_exact);
    std::cout << "b0=" << b[0] << " b1=" << b[1] << "\n";   // 期望 0.255782 / 3.47687

    auto x = thomas_solve(a, d, c, b);
    for (std::size_t i = 0; i < n; ++i)
        EXPECT_NEAR(x[i], x_exact[i], 1e-12);
}
