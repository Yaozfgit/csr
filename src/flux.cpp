#include <csr/flux.hpp>
#include <cmath>

namespace csr {
double rusanov_flux(double uL, double uR, double a) {
    const double fL = a * uL;
    const double fR = a * uR;
    const double s_max = std::abs(a);
    return 0.5* (fL + fR )-0.5*s_max*(uR-uL);
    }
}