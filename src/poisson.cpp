#include "csr/poisson.hpp"

CSRMatrix build_poisson(std::size_t n) {
    CSRMatrix A(n);
    double h2inv = static_cast<double>(n+1)*(n+1); //1/h =(n+1)
    for (std::size_t i = 0; i<n; ++i) {
        if (i >0) A.insert(i, i-1, -1.0 * h2inv);
        A.insert(i, i, 2.0 * h2inv);
        if (i+1 < n) A.insert(i, i+1,-1.0*h2inv);
    }
    return A;
}
