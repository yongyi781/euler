#pragma once

#include "../modular_arithmetic.hpp"
#include "../prime.hpp"
#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates the prime factorization of a number.
template <integral2 T, std::ranges::view V = std::ranges::empty_view<int>> class factor : public it_base
{
  public:
    using value_type = PrimePower<T>;

    factor() = default;
    constexpr explicit factor(T n, V spfs = {}) : _n(std::move(n)), _spfs(std::move(spfs))
    {
        assert(_n != 0 && "0 does not have a factorization");
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        T n = _n;
        if (n < 0)
        {
            if (!callbackResult(f, PrimePower<T>{-1, 1}))
                return result_break;
            n = -n;
        }
        T start = 2;
        while (n > 1)
        {
            T p = !_spfs.empty() && n < (int64_t)_spfs.size() ? _spfs[(size_t)n] : smallestPrimeFactor(n, start);
            if (!callbackResult(f, PrimePower<T>{p, valuationDivide<true>(n, p)}))
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
    V _spfs;
};

template <integral2 T, std::ranges::random_access_range Range>
factor(T, Range &&) -> factor<T, std::views::all_t<Range>>;
} // namespace it

/// Gives the prime factorization of a number.
template <integral2 T, std::ranges::range Range = std::vector<int>>
constexpr Factorization<T> factor(T num, Range &&spfs = {})
{
    Factorization<T> result;
    result.reserve(8);
    it::factor(std::move(num), std::forward<Range>(spfs))([&](auto &&pe) { result.push_back(pe); });
    return result;
}

/// Calculates Euler's totient function.
template <integral2 T> constexpr T totient(T n) // Pass by value intentional
{
    it::factor{n}([&](auto &&pe) { n -= n / pe.first; });
    return n;
}

/**
 * Calculates the multiplicative order of `a` modulo `modulus` with respect to the given prime factorization of the
 * totient.
 *
 * @param a The base number.
 * @param modulus The modulus.
 * @param totient The totient of `modulus`.
 * @return The multiplicative order of `a` modulo `modulus`.
 */
template <integral2 Ta, integral2 Tp, integral2 Tt, std::ranges::range Range = std::vector<int>>
constexpr Tt multiplicativeOrder(const Ta &a, const Tp &modulus, Tt totient, Range &&spfs = {})
{
    using euler::gcd;
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    it::factor(totient, std::forward<Range>(spfs))([&](auto &&pe) {
        auto [q, _] = pe;
        while (totient % q == 0 && powmSafe(a, totient / q, modulus) == 1)
            totient /= q;
    });
    return totient;
}

/**
 * Calculates the multiplicative order of `a` modulo `modulus`.
 *
 * @param a The base number.
 * @param modulus The modulus.
 * @return The multiplicative order of `a` modulo `modulus`.
 */
template <integral2 Ta, integral2 Tp> constexpr Tp multiplicativeOrder(const Ta &a, const Tp &modulus)
{
    using euler::gcd;
    using std::gcd;
    if (gcd(a, modulus) != 1)
        return 0;
    return multiplicativeOrder(a, modulus, totient(modulus));
}
} // namespace euler
