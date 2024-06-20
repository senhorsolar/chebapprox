#include "ChebyshevN.h"
#include <iostream>

double f (double x, double y, double z)
{
    return x * y * z;
}

int main()
{
    std::array<double, 3> a = {-1.0, -1.0, -1.0};
    std::array<double, 3> b = {1.0, 1.0, 1.0};
    std::size_t N = 10;
    std::size_t degree = 5;
    auto func = std::function (f);
    Chebyshev::ChebyshevN cheb (func, a, b); // can exclude a and b
    cheb.Fit (degree, N);

    std::array<double, 3> fargs = {0.5, 0.5, 0.5};

    std::cout << "Actual: " << func (0.5, 0.5, 0.5) << '\n';
    std::cout << "Approx: " << cheb.Approximate (0.5, 0.5, 0.5) << '\n';
    std::cout << "Approx alt call: " << cheb.Approximate (fargs) << '\n';
}
