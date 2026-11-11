#include <csr/poisson.hpp>
#include <csr/grid.hpp>

CSRMatrix build_poisson(const Grid1D& g) {    
    const std::size_t n = g.n_interior();    
    const double h2inv = 1.0 / (g.h() * g.h());   
    CSRMatrix A(n); 
    for (std::size_t i = 0; i<n; ++i) {
        const double diag = (i == 0 || i == n-1) ? 3.0 * h2inv
                                                 : 2.0 * h2inv;
        if (i >0) A.insert(i, i-1, -h2inv);
        A.insert(i, i, diag);
        if (i+1 < n) A.insert(i, i+1, -h2inv);
    }
    return A;
}
