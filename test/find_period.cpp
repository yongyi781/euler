#include "euler/ZMod.hpp"
#include "euler/algorithm.hpp"

#include "common.hpp"

using namespace std;
using Int = int64_t;
using Z = ZMod<(Int)1'000'000'007>;

/// Shortcut to write one-line lambdas
#define fun(var, ...) ([&](auto &&var) -> decltype(auto) { return __VA_ARGS__; })

int main()
{
    Int const preperiod = 39911;
    Int const period = 21353;
    auto periodResult = findPeriod(fun(a, a * a + 2), Z(0));
    assertEqual(periodResult.period, period);
    assertEqual(periodResult.preperiod, preperiod);
    Int const N = preperiod + 2 * period;
    auto v = unfold(N + 1, Z(0), fun(a, a * a + 2));
    auto fn = [&](Int n) { return v[n]; };
    for (Int n = preperiod + 1; n <= preperiod + 2 * period; ++n)
    {
        // Z a = 0;
        auto res1 = sumPeriodic([&](Int n) { return v[n]; }, preperiod, period, n / 2, n);
        auto res2 = sum(n / 2, n, fn);
        assertEqual(res1, res2);
    }
    pass("findPeriod");
}
