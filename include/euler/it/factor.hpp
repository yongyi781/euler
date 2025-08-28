#pragma once

#include "../PF.hpp"
#include "../modular_arithmetic.hpp"
#include "../prime.hpp"
#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates the prime factorization of a number.
template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>> class factor : public it_base
{
    T _n;
    std::reference_wrapper<const SPFSieve> _spfs;

  public:
    using value_type = PrimePower<T>;

    factor() = default;
    constexpr factor(T n, const SPFSieve &spfs = {}) : _n(std::move(n)), _spfs(std::ref(spfs))
    {
        assert(_n != 0 && "0 does not have a factorization");
    }

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        T n = _n;
        if (n < 0)
        {
            if (!callbackResult(f, PrimePower<T>{T(-1), T(1)}))
                return result_break;
            n = -n;
        }
        T start = 2;
        while (n > 1)
        {
            T p = !_spfs.get().empty() && n < (int64_t)_spfs.get().size() ? _spfs.get()[(size_t)n]
                                                                          : smallestPrimeFactor(n, start);
            if (!callbackResult(f, PrimePower<T>{p, removeFactors<true>(n, p)}))
                return result_break;
            if (p == 2)
                start = 3;
            else
                start = p + 2;
        }
        return result_continue;
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
template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>>
constexpr PF<T> factor(T num, SPFSieve &&spfs = {})
{
    PF<T> result;
    result.data().reserve(8);
    it::factor(std::move(num), spfs).appendTo(result.data());
    return result;
}

/// Calculates Euler's totient function.
template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>>
constexpr T totient(T n, SPFSieve &&spfs = {}) // Pass by value intentional
{
    it::factor{n, spfs}([&](auto &&pe) { n -= n / pe.first; });
    return n;
}

/// Calculates the multiplicative order of `a` modulo `modulus` with respect to the given prime factorization of the
/// totient.
///
/// @param a The base number.
/// @param modulus The modulus.
/// @param totient The totient of `modulus`.
/// @return The multiplicative order of `a` modulo `modulus`.
template <integral2 Ta, integral2 Tp, integral2 Tt, typename SPFSieve = std::vector<int>>
constexpr Tt multiplicativeOrder(const Ta &a, const Tp &modulus, Tt totient, SPFSieve &&spfs = {})
{
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    it::factor(totient, std::forward<SPFSieve>(spfs))([&](auto &&pe) {
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
template <integral2 Ta, integral2 Tp> constexpr Tp multiplicativeOrder(const Ta &a, const Tp &modulus)
{
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    return multiplicativeOrder(a, modulus, totient(modulus));
}

/// Returns (s, t) where d = s²t and t is squarefree.
template <integral2 T, std::ranges::range Range> constexpr std::pair<T, T> sqfreeDecompose(T num, Range &&factorization)
{
    T s = 1;
    for (auto &&[p, e] : std::forward<Range>(factorization))
    {
        T const q = pow(p, e / 2);
        s *= q;
        num /= q * q;
    }
    return {s, num};
}

/// Returns (s, t) where d = s²t and t is squarefree.
template <integral2 T> constexpr std::pair<T, T> sqfreeDecompose(T num)
{
    T s = 1;
    it::factor{num}([&](auto &&pe) {
        auto &&[p, e] = pe;
        T const q = pow(p, e / 2);
        s *= q;
        num /= q * q;
    });
    return {s, num};
}

/// Returns the radical of the given number.
template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>> T radical(T num, SPFSieve &&spfs = {})
{
    return it::factor{num, spfs}.map([](auto &&pe) { return pe.first; }).product();
}

/// Function to find smallest primitive root of p.
template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>>
constexpr T primitiveRoot(T p, SPFSieve &&spfs = {})
{
    if (p == 2)
        return 1;
    T const phi = p - 1;
    auto const pf = factor(phi, std::forward<SPFSieve>(spfs));
    for (T r = 2; r <= phi; ++r)
    {
        bool found = true;
        for (const auto &[q, e] : pf)
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

/// Returns the sum of a function at integers coprime to the given integer in the range [1, limit].
template <typename SummatoryFun, integral2 Tk, typename T> constexpr auto sumCoprime(SummatoryFun F, Tk k, T limit)
{
    thread_local std::vector<Tk> primes;
    primes.clear();
    it::factor(k).map([&](auto &&t) { return t.first; }).appendTo(primes);
    return sumCoprime(std::move(F), primes.begin(), primes.end(), limit);
}

/// Returns the number of integers coprime to the given integer in the range [1, limit].
template <integral2 Tk, typename T> constexpr T countCoprime(Tk k, T limit)
{
    thread_local std::vector<Tk> primes;
    primes.clear();
    it::factor(k).map([&](auto &&t) { return t.first; }).appendTo(primes);
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
template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>>
constexpr std::vector<T> sortedDivisors(T n, SPFSieve &&spfs = {})
{
    std::vector<T> res{T(1)};
    it::factor(n, spfs)([&](auto &&pe) { sortedDivisorsMerge(res, pe.first, pe.second); });
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
