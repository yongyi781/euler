#pragma once

#include "../prime.hpp"
#include "factor.hpp"

inline namespace euler
{
namespace it
{
/// Enumerate the divisors of a number.
///
/// Usage: `divisors(factorization)`. If the number `n` hasn't been factored yet, you probably want to type
/// `divisors(factor(n))`.
template <std::ranges::view V> class divisors_t : public it_base
{
  public:
    using value_type = std::ranges::range_value_t<V>::first_type;

    divisors_t() = default;
    constexpr divisors_t(V factorization) : _factorization(std::move(factorization)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _enumerate(std::ranges::begin(_factorization), std::ranges::end(_factorization), value_type(1),
                          std::move(f));
    }

    /// Returns the number of divisors, using a faster algorithm. Make sure to use a larger type `T`
    /// if you expect an answer larger than `2^64 - 1`.
    template <typename T = std::size_t> [[nodiscard]] constexpr T size() const
    {
        return it::wrap(std::ranges::ref_view(_factorization))
            .map([](auto &&pe) -> T { return pe.second + 1; })
            .product();
    }

    /// Returns sum of the divisors, using a faster algorithm.
    template <typename T = value_type> [[nodiscard]] constexpr T sum() const
    {
        return it::wrap(std::ranges::ref_view(_factorization))
            .map([](auto &&pe) -> T {
                if (pe.second == 1)
                    return pe.first + 1;
                return T(pow(T(pe.first), pe.second + 1) - 1) / T(pe.first - 1);
            })
            .product();
    }

  private:
    using It = std::ranges::iterator_t<const V>;
    V _factorization;

    template <std::invocable<value_type> Fun>
    constexpr result_t _enumerate(It it, It end, const value_type &current, Fun f) const
    {
        if (!callbackResult(f, current))
            return result_break;
        for (; it != end; ++it)
        {
            auto &&[p, e] = *it;
            value_type n = current;
            for (int i = 1; i <= e; ++i)
            {
                n *= p;
                if (!_enumerate(std::ranges::next(it), end, n, f))
                    return result_break;
            }
        }
        return result_continue;
    }
};

template <std::ranges::view V> it::divisors_t<V> divisors(V factorization) { return {std::move(factorization)}; }

template <std::ranges::range Range> it::divisors_t<std::views::all_t<Range>> divisors(Range &&r)
{
    return {std::forward<Range>(r)};
}

template <integral2 T> it::divisors_t<std::views::all_t<Factorization<T>>> divisors(T num)
{
    return {factor(std::move(num)).to()};
}
} // namespace it

/// Generates a vector of divisors based on the given factorization.
template <std::ranges::range Range> constexpr auto divisors(Range &&factorization)
{
    using T = std::ranges::range_value_t<Range>::first_type;
    auto d = it::divisors(std::forward<Range>(factorization));
    std::vector<T> result(d.size());
    auto it = result.begin();
    d([&](auto &&d) { *it++ = std::forward<decltype(d)>(d); });
    return result;
}

template <integral2 T, typename SPFSieve = std::ranges::empty_view<T>>
constexpr std::vector<T> divisors(T num, SPFSieve &&spfs = {})
{
    if (num == 1)
        return {1};
    return divisors(factor(std::move(num), std::forward<SPFSieve>(spfs)));
}
} // namespace euler
