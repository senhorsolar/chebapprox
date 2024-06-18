/*
 * Chebyshev1 Implementation
 */

#include "Chebyshev1.h"

namespace Chebyshev {

Chebyshev1::Chebyshev1 (std::size_t N)
    : m_N (N),
      m_fitted (false),
      m_a (0),
      m_b (0),
      m_coeffs (N),
      m_f ()
{
}

std::vector<double>
Chebyshev1::ChebyshevPolynomials (double x)
{
    std::vector<double> T (m_N);
    if (m_N > 0)
        T[0] = 1;
    if (m_N > 1)
        T[1] = x;
    if (m_N > 2) {
        for (std::size_t i = 2; i < m_N; ++i) {
            T[i] = 2*x*T[i-1] - T[i-2];
        }
    }
    return T;
}

void
Chebyshev1::Fit (std::function<double(double)> f, double a, double b)
{
    m_f = f;
    m_a = a;
    m_b = b;

    std::vector<std::vector<double>> all_polys;
    std::vector<double> y;

    for (std::size_t k = 1; k <= m_N; ++k) {
        double u = std::cos ((2*k - 1)/(2.0*m_N) * std::numbers::pi);
        double x = (b-a)/2.0 * u + (a+b)/2.0;
        all_polys.push_back (this->ChebyshevPolynomials (u));
        y.push_back (f (x));
    }

    for (std::size_t j = 0; j < m_N; ++j) {
        double coeff = 0.0;
        for (std::size_t i = 0; i < m_N; ++i) {
            coeff += y[i] * all_polys[i][j];
        }
        m_coeffs[j] = 2.0/m_N * coeff;
    }

    m_fitted = true;
}

double
Chebyshev1::Approximate (double x)
{
    if (!m_fitted) {
        std::cerr << "Chebyshev1 not fitted yet, returning 0.\n";
        return 0.0;
    }
    if (m_N < 1) {
        return 0.0;
    }

    double u = (2*x - m_a - m_b) / (m_b - m_a);
    std::vector<double> polys = this->ChebyshevPolynomials (u);

    return std::inner_product(
        m_coeffs.begin (), m_coeffs.end (), polys.begin (), -0.5*m_coeffs[0]);
}

}
