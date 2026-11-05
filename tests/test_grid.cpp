#include <csr/grid.hpp>
#include <gtest/gtest.h>

TEST(Grid1DTest, Spacing) {
    Grid1D g(9);  //内部9个点
    EXPECT_EQ(g.n_interior(),9u);
    EXPECT_EQ(g.n_total(),11u);
    EXPECT_NEAR(g.h(),0.1,1e-15);
}
TEST(Grid1DTest, Endpoints) {               // G2：边界点
    Grid1D g(3);
    EXPECT_NEAR(g.x(0), 0.0, 1e-15);
    EXPECT_NEAR(g.x(4), 1.0, 1e-15);        // n+1 = 4
}
TEST(Grid1DTest, InteriorConsistent) {      // G3：内部点 = (i)/(n+1)
    Grid1D g(9);
    for (std::size_t i = 1; i <= 9; ++i)
        EXPECT_NEAR(g.x(i), i * 0.1, 1e-12);
}
TEST(Grid1DTest, NoDuplicateSources) {      // G4：与收敛测试用同一来源
    Grid1D g(99);
    EXPECT_NEAR(g.interior_points()[0], 0.01, 1e-12);
}