#include "Chebyshev1.h"

double f (double x)
{
    return std::cos (x) + 0.3*std::pow (x, 3) + 2*std::pow (x, 2) + x - 10;
}

int main()
{
    std::size_t N = 7;
    Chebyshev::Chebyshev1 cheb (N);

    cheb.Fit (f);

    double x = 0.5;
    double z = f (0.5);
    double z_approx = cheb.Approximate (x);

    std::cout << "f(x): " << z << ", approx(x): " << z_approx << '\n';
}
