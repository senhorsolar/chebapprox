#include "ChebyshevN.h"
#include <iostream>
#include <tuple>
#include <type_traits>

// template<double... D>
// void f(D... args)
// {

// }


template<typename T>
concept Float = std::is_floating_point_v<T>;

template<Float R, Float... Args>
R func (Args... args)
{
    constexpr std::size_t nargs = sizeof...(args);
    std::cout << "Called with " << nargs << " args\n";
    for (const auto arg : {args...}) {
        std::cout << "--arg : " << arg << '\n';
    }
    return 0.0;
}

// Foo<double> f1;
// Foo<double, double> f2;
// Foo<double, double, int> f3;

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
    auto f2 = std::function (f);
    Chebyshev::ChebyshevN cheb (f2);//, a, b);
    cheb.Fit (degree, N);

    std::array<double, 3> fargs = {0.5, 0.5, 0.5};

    std::cout << "Actual: " << f2 (0.5, 0.5, 0.5) << '\n';
    std::cout << "Approx: " << cheb.Approximate (0.5, 0.5, 0.5) << '\n';
    std::cout << "V Approx: " << cheb.ApproximateV (fargs) << '\n';
}
