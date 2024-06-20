/*
 * Chebyshev1 Header
 */

#ifndef CHEBYSHEV_1_H_
#define CHEBYSHEV_1_H_

#include <cmath>
#include <functional>
#include <numbers>
#include <numeric>
#include <vector>
#include <iostream>

namespace Chebyshev {

/*
 * \class
 *
 * Class for approximating univariate functions
 */
class Chebyshev1
{
    public:
        /*
         * \brief ctor
         * \param[in] f Univariate function to approximate
         * \param[in] a Lower bound of domain(f) (default -1)
         * \param[in] b Upper bound of domain(f) (default 1)
         */
        Chebyshev1 (std::function<double(double)> f,
                    double a = -1.0, double b = 1.0);

        /*
         * \brief Find polynomial coefficients to best approximate function f
         * \param degree Polynomial degree
         */
        void Fit (std::size_t degree);

        /*
         * \brief Approximate f(x)
         * \param[in] x Input value
         * \return ~f(x)
         */
        double Approximate (double x);

        /*
         * \brief Calculate chebyshev polynomials at input x
         * \param[in] x Input value
         * \return T[i](x) for i=0...N
         */
        std::vector<double> ChebyshevPolynomials (double x, std::size_t degree);

    private:
        std::function<double(double)> m_f; ///< Function
        double m_a; ///< Lower bound of domain(f)
        double m_b; ///< Upper bound of domain(f)
        std::size_t m_degree; ///< Polynomial degree
        bool m_fitted; ///< flag to indicate if fitted or not
        std::vector<double> m_coeffs; ///< Fitted poly coefficients
};

} // namespace

#endif // CHEBYSHEV_1_H_
