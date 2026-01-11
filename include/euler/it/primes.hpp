#pragma once

#include "base.hpp"
#include "euler/PF.hpp"
#include <primesieve/iterator.hpp>

namespace euler
{
namespace it
{
/// Enumerates primes, using the primesieve library.
class primes : public it_base
{
    u64 start_ = 2;
    u64 stop_ = std::numeric_limits<u64>::max();

  public:
    using value_type = u64;

    primes() = default;
    constexpr primes(u64 start, u64 stop = std::numeric_limits<u64>::max()) : start_(start), stop_(stop) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if (start_ > stop_)
            return result_continue;
        primesieve::iterator it(start_, stop_);
        for (u64 p = it.next_prime(); p <= stop_; p = it.next_prime())
            if (!callbackResult(f, p))
                return result_break;
        return result_continue;
    }
};
} // namespace it

/// Creates a boolean vector of 0 .. `limit` with prime indices set to true.
inline std::vector<bool> primeSieve(size_t limit)
{
    std::vector<bool> res(limit + 1);
    it::primes(2, limit)([&](auto &&p) { res[p] = true; });
    return res;
}

/// Generates a sieve of squarefree numbers up to a given limit.
template <typename T = bool> std::vector<T> squarefreeSieve(size_t limit)
{
    std::vector<T> sieve(limit + 1, true);
    sieve[0] = false;
    size_t const s = isqrt(limit);
    it::primes(2, s)([&](u64 p) {
        for (size_t j = p * p; j <= limit; j += p * p)
            sieve[j] = false;
    });
    return sieve;
}

/// Linear non-parallel sieve for the totient function. O(n).
template <std::integral T> constexpr std::vector<T> totientSieve_seq(T limit)
{
    std::vector<T> phi(limit + 1);
    std::vector<T> primes;
    phi[1] = 1;
    primes.reserve(limit / std::max(1.0, log(limit)));

    for (T i = 2; i <= limit; i++)
    {
        if (phi[i] == 0)
        {
            phi[i] = i - 1;
            primes.push_back(i);
        }
        for (T p : primes)
        {
            if (!mulLeq(i, p, limit))
                break;
            if (i % p == 0)
            {
                phi[i * p] = phi[i] * p;
                break;
            }
            phi[i * p] = phi[i] * (p - 1);
        }
    }
    return phi;
}

/// Factors `n!`.
template <integral2 T = int> constexpr PF<T> factorFactorial(int n)
{
    return it::primes(2, n).map([&](auto &&p) { return PrimePower<T>{p, factorialValuation(n, p)}; }).to();
}

/// Factors `binomial(n, r)`.
template <integral2 T = int> constexpr PF<T> factorBinomial(int n, int r)
{
    return it::primes(2, n)
        .map([&](auto &&p) {
            auto e = factorialValuation(n, p) - factorialValuation(n - r, p) - factorialValuation(r, p);
            return PrimePower<T>{p, e};
        })
        .to();
}
} // namespace euler
