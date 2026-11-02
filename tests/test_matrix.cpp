#include <gtest/gtest.h>
#include <csr/matrix.hpp>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

//----------测试1：矩阵构造与元素访问
TEST(MatrixTest, ConstructAndAccess) {
    Matrix m(3);
    EXPECT_EQ(m.size(), 3u);
    m(1, 1) = 5.0;
    EXPECT_DOUBLE_EQ(m(1, 1), 5.0);
    EXPECT_THROW(m(3,0), std::out_of_range);  //越界抛出异常
}

//---------测试2：SpMV vs 手工计算
TEST(MatrixTest, SpMVSimple) {
    Matrix m(2);
    m(0,0) = 1; m(0,1) = 2;
    m(1,0) = 3; m(1,1) = 4; 
    auto y = m.multiply(std::vector<double>{1.0,1.0});
    EXPECT_DOUBLE_EQ(y[0],3.0);
    EXPECT_DOUBLE_EQ(y[1],7.0);
}

//---------测试3：SpMV 对基向量 e_i = 矩阵的第i列
TEST(MatrixTest, SpMVUnitVectors) {
    Matrix m(3);
    m(0,0) = 1; m(0,1) = 2; m(0,2) = 3;
    m(1,0) = 4; m(1,1) = 5; m(1,2) = 6;
    m(2,0) = 7; m(2,1) = 8; m(2,2) = 9;
    for (std::size_t j = 0; j<3; ++j) {
        std::vector<double> e(3,0.0);
        e[j] = 1.0;
        auto y = m.multiply(e);
        EXPECT_DOUBLE_EQ(y[0],m(0,j));
        EXPECT_DOUBLE_EQ(y[1],m(1,j));
        EXPECT_DOUBLE_EQ(y[2],m(2,j));
    }
}

//--------测试4： 物理问题——二阶袋鼠（差分vs 解析）
TEST(MatrixTest, SecondDerivativeLaplacian) {
    const std::size_t n = 64;
    const double k = 2.0;
    const double h = M_PI /(n-1);
    Matrix L(n);
    for (std::size_t i = 0; i< n; ++i) {
        L(i,i) = -2.0;
        if (i>0) L(i, i-1) = 1.0;
        if (i< n-1) L(i,i+1) = 1.0;
    }
    std::vector<double> f(n);
    for (std::size_t i=0; i < n; ++i) f[i] = std::sin(k*i*h);

    auto d2 = L.multiply(f);

    //差分阶段误差O(h^2),容差1e-2是理论合理值
    for (std::size_t i = 1; i< n-1; ++i) {
        double exact = -k * k * std::sin(k*i*h);
        EXPECT_NEAR(d2[i] / (h * h), exact, 1e-2);
    }
}

//----------测试5:1x1 边界情况--------
TEST(MatrixTest, OneByOne) {
    Matrix one(1);
    one(0,0) = 3.0;
    auto y = one.multiply(std::vector<double>{2.0});
    EXPECT_DOUBLE_EQ(y[0], 6.0);    
}

//----------测试6:全负矩阵的max--------
TEST(MatrixTest, MaxAllNegative) {
    Matrix neg(2);
    neg(0,0) = -5.0;  neg(0,1) = -2.0;
    neg(1,0) = -8.0; neg(1,1) = -3.0;
    EXPECT_DOUBLE_EQ(neg.max(), -2.0);
}

/*
int main() {
    //测试1：单位矩阵 * x = x
    Matrix I(4);
    for (std::size_t i = 0; i < 4; ++i) I (i, i) = 1.0;  //4*4的单位矩阵（每个元素都是1）
    std::vector<double> x = {1.0,2.0,3.0,4.0};
    auto y = I.multiply(x);
    for (int i = 0; i < 4; ++i) assert (std::fabs(y[i] - x [i])< 1e-14);

    //测试2：越界了要抛出异常
    bool threw = false;
    try {I(4,0) =1.0; } //n=4时，合法下标是0~3
    catch (const std::out_of_range&) {threw = true; }
    assert (threw);

    //测试3：维度不匹配要抛出异常
    threw = false;
    std::vector<double> bad = {1.0,2.0};
    try { I.multiply(bad); }
    catch (const std::invalid_argument&) { threw = true; }
    assert (threw);
     
    //测试4（物理验证）：1D 二阶导数矩阵作用在 sin(kx) 上
    //二阶中心差分： f''(x) ≈ (f(x+h)-2f(x)+f(x-h))/h²
    //解析解：(sin(kx))'' = -k² sin(kx)
    const std:: size_t n = 2048;
    const double L = M_PI, h = L/(n-1), k = 2.0;
    Matrix D(n);
    for (std::size_t i = 1; i < n - 1; ++i) {
        D(i, i - 1) = 1.0; D(i,i)=-2.0;  D(i, i+1) = 1.0;
    }
    std::vector<double> f(n);
    for(std::size_t  i =0; i<n; ++i)  f[i] =std::sin(k * i * h);

    auto d2 = D.multiply(f);
    for (std::size_t i = 1; i< n-1; ++i) {
        double exact = -k * k *std:: sin(k * i * h);
        // 在 for 循环里临时加：
        //double err = std::fabs(d2[i] / (h * h) - exact);
       // if (err > 1e-4)
       // std::printf("i=%zu  numeric=%.6f  exact=%.6f  err=%.2e\n", i, d2[i] / (h * h), exact, err);

        assert (std::fabs(d2[i]/(h * h) - exact) < 1e-5);
    }

    std::printf("all 4 tests passed!\n");

    Matrix neg(2);
    neg(0,0) = -5.0; neg(0,1) =-2.0;
    neg(1,0) = -8.0; neg(1,1) =-3.0;
    assert(neg.max() == -2.0);
    std:: printf("test max passed!\n");

    Matrix one(1);
    one(0,0) = 3.0;
    auto y1 = one.multiply(std::vector<double>{2.0});
    assert(y1[0]== 6.0);
    
    return 0; 
}
*/  