#include <csr/burgers.hpp>
#include <csr/grid.hpp>
#include <csr/linalg.hpp>
#include <csr/matrix.hpp>
#include <csr/thomas.hpp>

BurgersSolver::BurgersSolver(const Grid1D& g, double nu)
    : g_(g), nu_(nu), n_(g.n_interior()) {}

CSRMatrix BurgersSolver::build_diffusion() const {
    const double h2inv = 1.0 / (g_.h() * g_.h());
    CSRMatrix L(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        const double diag = (i == 0 || i == n_-1) ? 3.0*h2inv : 2.0*h2inv;
        if (i > 0)     L.insert(i, i-1, -h2inv);
        L.insert(i, i, diag);
        if (i+1 < n_)  L.insert(i, i+1, -h2inv);
    }
    return L;
}

std::vector<double> BurgersSolver::convection_rhs(const std::vector<double>& u) const {
    const double h = g_.h();
    std::vector<double> rhs(n_, 0.0);
    // ghost：齐次边界 u=0 处 u_ghost = -u0（和 Poisson 首行同一套外推）
    auto face_flux = [&](double ul, double ur) {
        return (ul >= 0) ? 0.5*ul*ul : 0.5*ur*ur;   // 迎风
    };
    for (std::size_t i = 0; i < n_; ++i) {
        double ul = (i == 0)   ? 0.0 : u[i-1];  // 左格心（边界处通量为 u=0 → 0.5*0*0=0）
        double ur = (i == n_-1)? 0.0 : u[i+1];
        // 界面 i+1/2 由 u[i] 和 u[i+1] 决定；边界界面通量恒为 0（齐次 Dirichlet）
        double F_right = (i == n_-1) ? 0.0 : face_flux(u[i], u[i+1]);
        double F_left  = (i == 0)    ? 0.0 : face_flux(u[i-1], u[i]);
        rhs[i] = -(F_right - F_left) / h;
    }
    return rhs;
}

std::vector<double> BurgersSolver::step(const std::vector<double>& u, double dt) const {
    std::vector<double> rhs = convection_rhs(u);

    std::vector<double> b(n_);
    for (std::size_t i = 0; i < n_; ++i)
        b[i] = u[i] + dt * rhs[i];        // 对流部分已显式算完   

    const double r = nu_*dt/(g_.h()*g_.h());

    std::vector<double> lower(n_,-r), diag(n_,1.0+2.0*r), upper(n_,-r);

    const std::size_t n = u.size();

    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
        diag[i] = 1.0 + 3.0 * r;   // 边界行：复刻原 CSR 的边界离散
        } else {
            diag[i] = 1.0 + 2.0 * r;   // 内部行
        }
    }
    lower[0] = 0.0; upper[n -1] = 0.0;

    return thomas_solve(lower, diag, upper, b);      
}