#include "euler/prime.hpp"
#include <euler/Quadratic.hpp>

#include "common.hpp"

using namespace std;

inline auto testModpe(i64 a, i64 b, i64 c, i64 p)
{
    i64 const d = b * b - 4 * a * c;
    int const he = log(1e4) / log(p);
    // int const he = 2;
    for (int e = 1; e <= he; ++e)
    {
        i64 const n = pow(p, e);
        auto const expected = it::range(0, n - 1).filter([&](i64 x) { return (a * x * x + b * x + c) % n == 0; }).to();
        auto const actual = Quadratic(a, b, c).solveModPrimePower(p, e);
        assertEqual(actual, expected);
    }
}

inline auto testModn(i64 a, i64 b, i64 c, i64 n)
{
    i64 const d = b * b - 4 * a * c;
    auto const expected = it::range(0, n - 1).filter([&](i64 x) { return (a * x * x + b * x + c) % n == 0; }).to();
    auto const actual = Quadratic(a, b, c).solveMod(n);
    assertEqual(actual, expected);
}

int main()
{
    for (i64 const p : primeRange(2, 7))
        for (i64 a = -6; a <= 6; ++a)
            for (i64 b = -6; b <= 6; ++b)
                for (i64 c = -6; c <= 6; ++c)
                    if (gcd({a, b, c}) == 1)
                        testModpe(a, b, c, p);
    for (i64 n = 2; n <= 20; ++n)
        for (i64 a = -6; a <= 6; ++a)
            for (i64 b = -6; b <= 6; ++b)
                for (i64 c = -6; c <= 6; ++c)
                    if (gcd({a, b, c}) == 1)
                        testModn(a, b, c, n);
}
