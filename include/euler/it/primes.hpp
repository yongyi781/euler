#pragma once

#include "base.hpp"
#include "euler/PF.hpp"
#include <primesieve/iterator.hpp>

inline namespace euler
{
namespace it
{
/// Enumerates primes, using the primesieve library.
class primes : public it_base
{
  public:
    using value_type = int64_t;

    primes() = default;
    constexpr primes(int64_t start, int64_t stop = std::numeric_limits<int64_t>::max()) : _start(start), _stop(stop) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if (_start > _stop)
            return result_continue;
        primesieve::iterator it(
            _start, _stop == std::numeric_limits<int64_t>::max() ? std::numeric_limits<uint64_t>::max() : _stop);
        for (int64_t p = it.next_prime(); p <= _stop; p = it.next_prime())
            if (!callbackResult(f, p))
                return result_break;
        return result_continue;
    }

  private:
    int64_t _start = 2;
    int64_t _stop = std::numeric_limits<int64_t>::max();
};
} // namespace it

/// Creates a boolean vector of 0 .. `limit` with prime indices set to true.
inline std::vector<bool> primeSieve(size_t limit)
{
    std::vector<bool> res(limit + 1, false);
    it::primes(2, limit)([&](auto &&p) { res[p] = true; });
    return res;
}

/// Factors `n!`.
template <integral2 T = int> constexpr PF<T> factorFactorial(int n)
{
    return it::primes{2, n}.map([&](auto &&p) { return PrimePower<T>{p, factorialValuation(n, p)}; }).to();
}

/// Factors `binomial(n, r)`.
template <integral2 T = int> constexpr PF<T> factorBinomial(int n, int r)
{
    return it::primes{2, n}
        .map([&](auto &&p) {
            auto e = factorialValuation(n, p) - factorialValuation(n - r, p) - factorialValuation(r, p);
            return PrimePower<T>{p, e};
        })
        .to();
}
} // namespace euler
