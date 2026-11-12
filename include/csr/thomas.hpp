#pragma once
#include <vector>
#include <stdexcept>

// 求解三对角线性系统，不修改输入
// a: 下副对角 (a[0] 未用), d: 主对角, c: 上副对角 (c[n-1] 未用), b: 右端项
// 返回解向量
inline std::vector<double> thomas_solve(const std::vector<double>& a,
                                        const std::vector<double>& d,
                                        const std::vector<double>& c,
                                        const std::vector<double>& b) {
    const std::size_t n = d.size();  //n是主对角线d的数量
    if (a.size() != n || c.size() != n || b.size() != n)
        throw std::invalid_argument("thomas_solve: size mismatch");
    std::vector<double> cp(n, 0.0), bp(n, 0.0), x(n, 0.0);

    //前向消元
    cp[0] = c[0]/d[0];  //cp的p是prime：上标
    bp[0] = b[0]/d[0];
    for (std::size_t i = 1; i<n; ++i) {
        const double m = d[i] - a[i]*cp[i-1];
        if (m == 0) throw std::runtime_error("thomas_solve: zero pivot");
        cp[i] = c[i]/m;
        bp[i] = (b[i]-a[i]*bp[i-1])/m;
    }

    //回代
    x[n-1] = bp[n-1];
    for (std::size_t i = n-1; i-->0;)
        x[i] = bp[i] - cp[i]*x[i+1];

    return x;
}                                       