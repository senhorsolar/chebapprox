/*
 * Chebyshev1 Implementation
 */

#include "Chebyshev1.h"

namespace Chebyshev {

Chebyshev1::Chebyshev1 (std::function<double(double)> f,
                        double a, double b)
    : m_f {f},
      m_a (a),
      m_b (b),
      m_degree (0),
      m_fitted (false),
      m_coeffs ()
{
}

std::vector<double>
Chebyshev1::ChebyshevPolynomials (double x, std::size_t degree)
{
    std::vector<double> T (degree + 1);
    T[0] = 1;
    if (degree > 0)
        T[1] = x;
    for (std::size_t k=1; k < degree; ++k) {
        T[k+1] = 2 * x * T[k] - T[k-1];
    }
    return T;
}

void
Chebyshev1::Fit (std::size_t degree)
{
    std::vector<std::vector<double>> all_polys;
    std::vector<double> y;

    std::size_t N = degree + 1;
    for (std::size_t k = 1; k <= N; ++k) {
        double u = std::cos ((2*k - 1)/(2.0*N) * std::numbers::pi);
        double x = (m_b-m_a)/2.0 * u + (m_a+m_b)/2.0;
        all_polys.push_back (this->ChebyshevPolynomials (u, degree));
        y.push_back (m_f (x));
    }

    m_coeffs.resize (N);
    for (std::size_t j = 0; j < N; ++j) {
        double coeff = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            coeff += y[i] * all_polys[i][j];
        }
        m_coeffs[j] = 2.0/N * coeff;
    }
    m_degree = degree;
    m_fitted = true;
}

double
Chebyshev1::Approximate (double x)
{
    if (!m_fitted) {
        std::cerr << "Chebyshev1 not fitted yet, returning 0.\n";
        return 0.0;
    }

    double u = (2*x - m_a - m_b) / (m_b - m_a);
    std::vector<double> polys = this->ChebyshevPolynomials (u, m_degree);

    return std::inner_product(
        m_coeffs.begin (), m_coeffs.end (), polys.begin (), -0.5*m_coeffs[0]);
}

}
