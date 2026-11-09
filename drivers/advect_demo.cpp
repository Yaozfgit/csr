// 一维线性对流演示：方波平移一个周期
// 用法: ./advect_demo [n] [cfl]
#include <csr/grid.hpp>
#include <csr/advect1d.hpp>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {
    const std::size_t n = (argc > 1) ? std::stoul(argv[1]) : 200;
    const double cfl    = (argc > 2) ? std::stod(argv[2])  : 0.9;

    Grid1D g(n);
    std::vector<double> u0(g.n_interior(), 0.0);
    for (std::size_t i = n/5; i < 2*n/5; ++i) u0[i] = 1.0;   // 方波 [0.2,0.4]

    csr::LinearAdvection adv(g, 1.0, u0);

    const int nsteps = static_cast<int>(1.0 / (cfl * g.h()));  // 与C1相同的凑整技巧
    const double dt = 1.0 / nsteps;

    std::cout << "n=" << n << " cfl=" << cfl
              << " steps=" << nsteps << " dt=" << dt << "\n";
    for (int s = 1; s <= nsteps; ++s) {
        adv.step(dt);
        if (s % (nsteps / 5) == 0) {          // 5 个时刻打印轮廓
            std::cout << "t=" << std::fixed << std::setprecision(2) << dt * s << " : ";
            for (std::size_t i = 0; i < g.n_interior(); i += n / 20)
                std::cout << std::setprecision(3) << adv.solution()[i] << " ";
            std::cout << "\n";
        }
    }
}
