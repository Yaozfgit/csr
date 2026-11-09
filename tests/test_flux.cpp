#include <gtest/gtest.h>
#include <csr/flux.hpp>

//F1： 光滑区域——左右相等时通量退化为物理通量 f=a*u
TEST(FluxTest, SmoothRegion) {
    EXPECT_NEAR(csr::rusanov_flux(1.0, 1.0, 1.0), 1.0, 1e-15);
}

//F2: 间断处——耗散项把通量往中间拉（激波出的行为）
TEST(FluxTest, JumpReducesFlux) {
    double f = csr::rusanov_flux(1.0, 0.0, 1.0); //0.5*(1+0)-0.5*(0-1)=1.0?
    EXPECT_NEAR(f, 1.0, 1e-15);  //手算验证：a=1时，0.5*（fL+fR)=0.5, 0.5*（uR-uL)=-0.5,合计为1.0
}

// F3: 对称性——齐次性（奇对称），通量变号
TEST(FluxTest, SignHomogeneity) {
    EXPECT_NEAR(csr::rusanov_flux(-0.8, -0.3, 1.0),
                -csr::rusanov_flux(0.8, 0.3, 1.0), 1e-15);
}

// F3-1: 对称性——齐次性（奇对称），a的反对称
TEST(FluxTest, MirrorSymmetry) {
    EXPECT_NEAR(csr::rusanov_flux(0.3, 0.8, -1.0),
                -csr::rusanov_flux(0.8, 0.3, 1.0), 1e-15);
}
