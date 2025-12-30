#pragma once

#include "../PF.hpp"
#include "../modular_arithmetic.hpp"
#include "../prime.hpp"
#include "base.hpp"
#include "euler/SPF.hpp"

namespace euler
{
namespace it
{
/// Enumerates the prime factorization of a number.
template <integral2 T, std::integral S = u32> class factor : public it_base
{
    T n_;
    const SPF<S> *spf_ = nullptr;

    /// Requires `n > 0`.
    template <typename Fun> constexpr result_t factorWithoutSpf(T n, Fun f) const
    {
        for (T start = 2; n != 1;)
        {
            T p = smallestPrimeFactor(n, start);
            if (!callbackResult(f, PrimePower<T>{p, removeFactors<true>(n, p)}))
                return result_break;
            if (p == 2)
                start = 3;
            else
                start = p + 2;
        }
        return result_continue;
    }

    /// Requires `spf_ != nullptr` and `0 < n < spf_->size()`.
    template <typename Fun> constexpr result_t factorWithSpf(T n, Fun f) const
    {
        for (T start = 2; n != 1 && n >= spf_->size();)
        {
            T p = smallestPrimeFactor(n, start);
            if (!callbackResult(f, PrimePower<T>{p, removeFactors<true>(n, p)}))
                return result_break;
            if (p == 2)
                start = 3;
            else
                start = p + 2;
        }
        S x = (S)std::move(n);
        while (x != 1)
            if (!callbackResult(f, spf_->remove(x)))
                return result_break;
        return result_continue;
    }

  public:
    using value_type = PrimePower<T>;

    constexpr explicit factor(T n) : n_(std::move(n)) { assert(n_ != 0 && "0 does not have a factorization"); }
    constexpr factor(T n, const SPF<S> &spf) : n_(std::move(n)), spf_(&spf)
    {
        assert(n_ != 0 && "0 does not have a factorization");
    }

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        T n = n_;
        if (n < 0)
        {
            if (!callbackResult(f, PrimePower<T>{T(-1), T(1)}))
                return result_break;
            n = -n;
        }
        return spf_ && !spf_->empty() ? factorWithSpf(std::move(n), std::move(f))
                                      : factorWithoutSpf(std::move(n), std::move(f));
    }

    /// Counts the number of divisors from the factorization.
    template <integral2 Z = T> constexpr Z countDivisors() const
    {
        return map([](auto &&pe) -> Z { return pe.second + 1; }).product();
    }

    /// Sums the divisors from the factorization.
    template <integral2 Z = T> constexpr Z sumDivisors() const
    {
        return map([](auto &&pe) -> Z {
                   auto &&[p, e] = pe;
                   if (e == 1)
                       return 1 + p;
                   if (e == 2)
                       return 1 + p + Z(p) * p;
                   return (pow(Z(p), e + 1) - 1) / Z(p - 1);
               })
            .product();
    }
};
} // namespace it

/// Gives the prime factorization of a number.
template <integral2 T, std::integral S = u32> constexpr PF<T> factor(T n, const SPF<S> &spf = {})
{
    PF<T> result;
    result.data().reserve(8);
    it::factor{std::move(n), spf}.appendTo(result.data());
    return result;
}

/// Calculates Euler's totient function.
template <integral2 T, std::integral S = u32>
constexpr T totient(T n, const SPF<S> &spf = {}) // Pass by value intentional
{
    it::factor{std::move(n), spf}([&](auto &&pe) { n -= n / pe.first; });
    return n;
}

/// Calculates the multiplicative order of `a` modulo `modulus` with respect to the given prime factorization of the
/// totient.
///
/// @param a The base number.
/// @param modulus The modulus.
/// @param totient The totient of `modulus`.
/// @return The multiplicative order of `a` modulo `modulus`.
template <integral2 Ta, integral2 Tp, integral2 Tt, std::integral S = u32>
constexpr Tt multiplicativeOrder(const Ta &a, const Tp &modulus, Tt totient, const SPF<S> &spf = {})
{
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    it::factor{totient, spf}([&](auto &&pe) {
        auto [q, e] = pe;
        while (e-- > 0 && powmSafe(a, totient / q, modulus) == 1)
            totient /= q;
    });
    return totient;
}

/// Calculates the multiplicative order of `a` modulo `modulus`.
///
/// @param a The base number.
/// @param modulus The modulus.
/// @return The multiplicative order of `a` modulo `modulus`.
template <integral2 Ta, integral2 Tp, std::integral S = u32>
constexpr Tp multiplicativeOrder(const Ta &a, const Tp &modulus, const SPF<S> &spf = {})
{
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    return multiplicativeOrder(a, modulus, totient(modulus, spf), spf);
}

/// Returns (s, t) where d = s²t and t is squarefree.
template <integral2 T, std::integral S = u32> constexpr std::pair<T, T> sqfreeDecompose(T n, const SPF<S> &spf = {})
{
    T s = 1;
    it::factor{n, spf}([&](auto &&pe) {
        auto &&[p, e] = pe;
        T const q = pow(p, e / 2);
        s *= q;
        n /= q * q;
    });
    return {s, n};
}

/// Returns the radical of the given number.
template <integral2 T, std::integral S = u32> T radical(T n, const SPF<S> &spf = {})
{
    return it::factor{std::move(n), spf}.map([](auto &&pe) { return pe.first; }).product();
}

/// Function to find smallest primitive root of p.
template <integral2 T> constexpr T primitiveRoot(T p, const PF<T> &fac_phi_p)
{
    if (p == 2)
        return 1;
    T const phi = p - 1;
    for (T r = 2; r <= phi; ++r)
    {
        bool found = true;
        for (const auto &[q, e] : fac_phi_p)
        {
            if (powmSafe(r, phi / q, p) == 1)
            {
                found = false;
                break;
            }
        }

        if (found)
            return r;
    }
    assert(false && "primitiveRoot: Should not reach here (maybe p wasn't prime).");
}

/// Function to find smallest primitive root of p.
template <integral2 T, std::integral S = u32> constexpr T primitiveRoot(T p, const SPF<S> &spf = {})
{
    if (p == 2)
        return 1;
    return primitiveRoot(p, factor(p - 1, spf));
}

/// Returns the sum of a function at integers coprime to the given integer in the range [1, limit].
template <typename Fun, typename SummatoryFun, integral2 Tk, typename T, std::integral S = u32>
constexpr auto sumCoprime(Fun f, SummatoryFun F, Tk k, T limit, const SPF<S> &spf = {})
{
    thread_local std::vector<Tk> primes;
    primes.clear();
    it::factor(k, spf).map([&](auto &&t) { return t.first; }).appendTo(primes);
    return sumCoprime(std::move(f), std::move(F), primes.begin(), primes.end(), limit);
}

/// Returns the number of integers coprime to the given integer in the range [1, limit].
template <integral2 Tk, typename T, std::integral S = u32>
constexpr T countCoprime(Tk k, T limit, const SPF<S> &spf = {})
{
    thread_local std::vector<Tk> primes;
    primes.clear();
    it::factor(k, spf).map([&](auto &&t) { return t.first; }).appendTo(primes);
    return countCoprime(primes.begin(), primes.end(), limit);
}

/// Performs one merge step for enumerating sorted divisors.
template <typename T, typename U> constexpr void sortedDivisorsMerge(std::vector<T> &res, U p, int e)
{
    auto tmp = res;
    while (e-- > 0)
    {
        for (auto &x : tmp)
            x *= p;
        res.append_range(tmp);
        std::ranges::inplace_merge(res, res.end() - tmp.size());
    }
}

/// Performs one merge step for enumerating sorted divisors, using a temporary vector passed in by the caller.
template <typename T, typename U>
constexpr void sortedDivisorsMerge(std::vector<T> &res, std::vector<T> &tmp, U p, int e)
{
    tmp = res;
    while (e-- > 0)
    {
        for (auto &x : tmp)
            x *= p;
        res.append_range(tmp);
        std::ranges::inplace_merge(res, res.end() - tmp.size());
    }
}

/// Returns the sorted list of divisors of `n`.
template <integral2 T, std::integral S = u32> constexpr std::vector<T> sortedDivisors(T n, const SPF<S> &spf = {})
{
    std::vector<T> res{T(1)};
    it::factor{n, spf}([&](auto &&pe) { sortedDivisorsMerge(res, pe.first, pe.second); });
    return res;
}

/// Returns the sorted list of divisors of `n` from its factorization.
template <std::ranges::range Range> constexpr auto sortedDivisors(Range &&factorization)
{
    using T = std::tuple_element_t<0, std::ranges::range_value_t<Range>>;
    std::vector<T> res{T(1)};
    for (auto &&[p, e] : factorization)
        sortedDivisorsMerge(res, p, e);
    return res;
}
} // namespace euler
