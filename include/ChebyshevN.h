/*
 * ChebyshevN Header
 */

#ifndef CHEBYSHEV_N_H_
#define CHEBYSHEV_N_H_

#include <array>
#include <cmath>
#include <functional>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct> // for tensor op

namespace Chebyshev {

template<typename T>
concept Float = std::is_floating_point_v<T>;

template<typename T, std::size_t N>
constexpr std::array<T, N> constArray (T c) {
    std::array<T, N> arr;
    arr.fill (c);
    return arr;
}

/*
 * \class
 *
 * Class for approximating multivariate functions
 */

template <Float R, Float... Args>
class ChebyshevN
{
    public:

        using Matrix=Eigen::Matrix<R, Eigen::Dynamic, Eigen::Dynamic>;
        using Vector=Eigen::Matrix<R, Eigen::Dynamic, 1>;

        /*
         * \brief ctor
         *
         * \param[in] f Multivariate function to approximate
         * \param[in] a Lower bound of domain(f)
         * \param[in] b Upper bound of domain(f)
         */
        ChebyshevN (
            std::function<R(Args...)> f,
            const std::array<R, sizeof...(Args)>& a = constArray<R, sizeof...(Args)> (-1.0),
            const std::array<R, sizeof...(Args)>& b = constArray<R, sizeof...(Args)> (1.0))
            : m_f {f},
              m_a (a),
              m_b (b),
              m_degree (0),
              m_fitted (false),
              m_coeffs ()
        {}

        /*
         * \brief Find polynomial coefficients to best approximate function f
         *
         * \param[in] degree Chebyshev polynomial degree
         * \param[in] N Number of sampling points (defaults to degree+1)
         */
        void Fit (std::size_t degree, std::size_t N=0);

        /*
         * \brief Approximate f(args...)
         * \param[in] args Input values
         * \return ~f(args...)
         */
        R Approximate (Args... args);

        template <class Tuple>
        R ApproximateV (Tuple&& t) {
            return std::apply ([&](auto...args) {
                return this->Approximate (std::forward<decltype(args)> (args)...);
            }, t);
        }

    protected:
        /*
         * \brief Calculate chebyshev polynomials at input x
         * \param[in] x Input value
         * \return T[i](x) for i=0...N
         */
        Matrix ChebyshevPolynomials (const Vector& x, std::size_t degree);
        Vector ChebyshevPolynomials (R x, std::size_t degree);
        static Vector ChebyshevNodes (std::size_t N);

    private:
        std::function<R(Args...)> m_f; ///< Function
        std::array<R, sizeof...(Args)> m_a;
        std::array<R, sizeof...(Args)> m_b;
        std::size_t m_degree; ///< max polynomial degree
        bool m_fitted; ///< flag to indicate if fitted or not
        Vector m_coeffs; ///< Fitted poly coefficients

};

template <Float R, Float... Args>
ChebyshevN<R, Args...>::Vector
ChebyshevN<R, Args...>::ChebyshevNodes (std::size_t N)
{
    Vector nodes = Vector::LinSpaced (N, 1, N).array () * 2 - 1;
    nodes.array () *= std::numbers::pi / (2.0*N);
    nodes.array () = nodes.array ().cos ();
    return nodes;
}

template <Float R, Float... Args>
ChebyshevN<R, Args...>::Matrix
ChebyshevN<R, Args...>::ChebyshevPolynomials (const Vector& x, std::size_t degree)
{
    Matrix T (degree+1, x.size ());
    T.row (0).array () = 1;
    if (degree > 0) {
        T.row (1).array () = x.transpose ().array ();
    }
    for (std::size_t k=1; k < degree; ++k) {
        T.row (k+1) =
            2 * x.transpose ().array () * T.row (k).array () - T.row (k-1).array ();
    }
    return T;
}

template <Float R, Float... Args>
ChebyshevN<R, Args...>::Vector
ChebyshevN<R, Args...>::ChebyshevPolynomials (R x, std::size_t degree)
{
    Vector xv (1);
    xv (0) = x;
    Vector poly = ChebyshevPolynomials (xv, degree);
    return poly;
}

template <Float R, Float... Args>
void ChebyshevN<R, Args...>::Fit (std::size_t degree, std::size_t N)
{
    if (N <= degree) {
        if (N > 0) {
            std::cout << "N <= degree, setting to degree+1\n";
        }
        N = degree + 1;
    }
    Matrix A;
    Vector nodes = this->ChebyshevNodes (N);
    Matrix poly = this->ChebyshevPolynomials (nodes, degree);
    for (std::size_t i=0; i < sizeof...(Args); ++i) {
        if (i == 0)
            A = poly;
        else {
            A = Eigen::KroneckerProduct (A, poly).eval ();
        }
    }
    //std::cout << "A dims: " << A.rows () << "," << A.cols () << '\n';

    std::array<Vector, sizeof...(Args)> arg_grid;
    for (std::size_t i=0; i < sizeof...(Args); ++i) {
        arg_grid[i] = Vector (N);
        for (std::size_t j=0; j < N; ++j) {
            arg_grid[i][j] = (m_b[i]-m_a[i])/2.0 * nodes[j] + (m_b[i]+m_a[i])/2.0;
            //arg_grid[i][j] = nodes[j];
        }
    }

    Vector z (A.cols ());
    for (std::size_t i {0};
         auto const& args : std::apply(std::views::cartesian_product, arg_grid)) {
        z (i++) = std::apply (m_f, args);
    }
    m_coeffs = A.transpose ().colPivHouseholderQr ().solve (z);
    //std::cout << "coeffs size: " << m_coeffs.size () << '\n';

    m_fitted = true;
    m_degree = degree;
}

template <Float R, Float... Args>
R ChebyshevN<R, Args...>::Approximate (Args... args)
{
    if (!m_fitted) {
        std::cerr << "Chebyshev polynomials not fitted yet, call Fit() first\n";
        return 0.0;
    }

    Vector A;
    for (std::size_t i=0; const auto arg : {args...}) {
        R u = (2 * arg - m_a[i] - m_b[i]) / (m_b[i] - m_a[i]);
        //R u = arg;
        if (i++ == 0)
            A = this->ChebyshevPolynomials (u, m_degree);
        else {
            A = Eigen::KroneckerProduct (
                A, this->ChebyshevPolynomials (u, m_degree)).eval ();
        }
    }

    // std::cout << "A size: " << A.rows () << "," << A.cols () << "\n";
    // std::cout << "Coefs size: " << m_coeffs.size () << '\n';
    return A.transpose ().dot (m_coeffs);
    //return 0.0;
}

} // namespace

#endif // CHEBYSHEV_N_H_
