#include "csr/matrix.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

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
    return 0; 
}