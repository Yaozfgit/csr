#include <csr/grid.hpp>
#include <gtest/gtest.h>

TEST(Grid1DTest, Spacing) {
    Grid1D g(9);  //内部9个点
    EXPECT_EQ(g.n_interior(),9u);
    EXPECT_EQ(g.n_total(),11u);
    EXPECT_NEAR(g.h(),1.0/9.0,1e-15);
    EXPECT_NEAR(g.h()*g.n_interior(),1.0,1e-15);
    EXPECT_DOUBLE_EQ(g.x(1), 0.5 * g.h());
}
TEST(Grid1DTest, GhostCells) {               // G2：边界点
    Grid1D g(9);
    EXPECT_NEAR(g.x(0), -0.5/9.0, 1e-15);
    EXPECT_NEAR(g.x(10), 1.0+0.5/9.0, 1e-15);        
}
TEST(Grid1DTest, CellCenters) {      // G3：内部点 = (i)/(n+1)
    Grid1D g(9);
    EXPECT_NEAR(g.x(1), 0.5/9.0, 1e-15);
    EXPECT_NEAR(g.x(9), 1.0-0.5/9.0, 1e-15);
    for (std::size_t i = 1; i <= 9; ++i)
        EXPECT_NEAR(g.x(i), (i-0.5) * g.h(), 1e-12);
}
TEST(Grid1DTest, NoDuplicateSources) {      // G4：与收敛测试用同一来源
    Grid1D g(99);
    EXPECT_NEAR(g.interior_points()[0], 0.5/99.0, 1e-12);
}