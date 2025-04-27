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
        return map([](auto &&pe) { return Z(pe.second + 1); }).product();
    }

    /// Sums the divisors from the factorization.
    template <integral2 Z = T> constexpr Z sumDivisors() const
    {
        return map([](auto &&pe) { return Z(pow(Z(pe.first), pe.second + 1) - 1) / Z(pe.first - 1); }).product();
    }

  private:
    T _n;
    std::reference_wrapper<const SPFSieve> _spfs;
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
    using euler::gcd;
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    it::factor(totient, std::forward<SPFSieve>(spfs))([&](auto &&pe) {
        auto [q, _] = pe;
        while (totient % q == 0 && powmSafe(a, totient / q, modulus) == 1)
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
    using euler::gcd;
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

/// Function to find smallest primitive root of p.
template <integral2 T, typename SPFSieve = std::vector<int>>
constexpr std::common_type_t<int64_t, T> primitiveRoot(const T &p, SPFSieve &&spfs = {})
{
    using Tp = std::common_type_t<int64_t, T>;
    if (p == 2)
        return Tp(1);
    T const phi = p - 1;
    auto const pf = factor(phi, std::forward<SPFSieve>(spfs));
    for (Tp r = 2; r <= phi; ++r)
    {
        bool flag = false;
        for (const auto &[q, e] : pf)
        {
            if (powm(r, phi / q, p) == 1)
            {
                flag = true;
                break;
            }
        }

        if (!flag)
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
} // namespace euler
