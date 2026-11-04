#include <csr/linalg.hpp>
#include <cmath>
#include <stdexcept>

std::vector<double> gauss_solve(Matrix A, std::vector<double> b) {
    const std::size_t n = b.size();
    for (std::size_t k = 0; k <n; ++k) {
        //选主元
        std::size_t p = k;
        for (std::size_t i = k+1; i<n; ++i)
            if (std::abs(A(i,k)) > std::abs(A(p,k))) p = i;
        if (std::abs(A(p,k)) < 1e-14)
            throw std::runtime_error("singular matrix");
        std::swap(b[k], b[p]);
        for (std::size_t j= k; j<n; ++j) std::swap(A(k, j),A(p,j));
        //消元
        for (std::size_t i=k+1; i<n; ++i) {
            double m = A(i,k)/A(k,k);
            for (std::size_t j=k; j<n; ++j) A(i,j) -= m*A(k,j);
            b[i] -= m *b[k];
        }
    }
    //回代
    std::vector<double> x(n);
    for (std::size_t i=n; i-- >0;) {
        double s = b[i];
        for (std::size_t j=i+1; j<n; ++j) s -= A(i,j)*x[j];
        x[i] = s/A(i,i);
    }
    return x;
}