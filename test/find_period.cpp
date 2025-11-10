#include "euler/ZMod.hpp"
#include "euler/algorithm.hpp"

#include "test.hpp"

using namespace std;
using namespace euler;
using Z = ZMod<1'000'000'007_i64>;

/// Shortcut to write one-line lambdas
#define fun(var, ...) ([&](auto &&var) -> decltype(auto) { return __VA_ARGS__; })

int main()
{
    i64 const preperiod = 39911;
    i64 const period = 21353;
    auto periodResult = findPeriod(fun(a, a * a + 2), Z(0));
    assertEqual(periodResult.period, period);
    assertEqual(periodResult.preperiod, preperiod);
    i64 const N = preperiod + 2 * period;
    auto v = unfold(N + 1, Z(0), fun(a, a * a + 2));
    auto fn = [&](i64 n) { return v[n]; };
    for (i64 n = preperiod + 1; n <= preperiod + 2 * period; ++n)
    {
        // Z a = 0;
        auto res1 = sumPeriodic([&](i64 n) { return v[n]; }, preperiod, period, n / 2, n);
        auto res2 = sum(n / 2, n, fn);
        assertEqual(res1, res2);
    }
}
