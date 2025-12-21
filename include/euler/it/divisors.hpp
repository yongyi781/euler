#pragma once

#include "../prime.hpp"
#include "factor.hpp"

namespace euler
{
namespace it
{
/// Enumerate the divisors of a number.
///
/// Usage: `divisors<T>(factorization)`. If the number `n` hasn't been factored yet, you probably want to type
/// `divisors<T>(factor(n))`.
template <integral2 T, std::ranges::view V> class divisors_t : public it_base
{
    using It = std::ranges::iterator_t<const V>;
    V _factorization;

    template <typename Fun> result_t enumerate(It it, It end, T cur, Fun f) const
    {
        if (!callbackResult(f, cur))
            return result_break;
        for (; it != end; ++it)
        {
            auto &&[p, e] = *it;
            T n = cur;
            for (int i = 1; i <= e; ++i)
            {
                n *= p;
                if (!enumerate(std::ranges::next(it), end, n, f))
                    return result_break;
            }
        }
        return result_continue;
    }

  public:
    using value_type = T;

    divisors_t() = default;
    constexpr divisors_t(V factorization) : _factorization(std::move(factorization)) {}

    template <typename Fun> result_t operator()(Fun f) const
    {
        return enumerate(std::ranges::begin(_factorization), std::ranges::end(_factorization), value_type(1),
                         std::move(f));
    }

    /// Returns the number of divisors, using a faster algorithm. Make sure to use a larger type `T`
    /// if you expect an answer larger than `2^64 - 1`.
    [[nodiscard]] constexpr T size() const
    {
        return it::wrap(std::ranges::ref_view(_factorization))
            .map([](auto &&pe) -> T { return pe.second + 1; })
            .product();
    }

    /// Returns sum of the divisors, using a faster algorithm.
    [[nodiscard]] constexpr T sum() const
    {
        return it::wrap(std::ranges::ref_view(_factorization))
            .map([](auto &&pe) -> T {
                if (pe.second == 1)
                    return pe.first + 1;
                return T(pow(T(pe.first), pe.second + 1) - 1) / T(pe.first - 1);
            })
            .product();
    }
};

template <integral2 T, std::ranges::range Range> it::divisors_t<T, std::views::all_t<Range>> divisors(Range &&fac)
{
    return {std::forward<Range>(fac)};
}

template <std::ranges::range Range>
it::divisors_t<typename std::ranges::range_value_t<Range>::first_type, std::views::all_t<Range>> divisors(Range &&fac)
{
    return {std::forward<Range>(fac)};
}

template <integral2 T> it::divisors_t<T, std::views::all_t<PF<T>>> divisors(T num)
{
    return {PF{factor(std::move(num)).to()}};
}
} // namespace it

/// Generates a vector of divisors based on the given factorization.
template <integral2 T, std::ranges::range Range> auto divisors(Range &&fac)
{
    auto d = it::divisors(std::forward<Range>(fac));
    std::vector<T> res((size_t)d.size());
    auto it = res.begin();
    d([&](auto &&d) { *it++ = std::forward<decltype(d)>(d); });
    return res;
}

/// Generates a vector of divisors based on the given factorization.
template <std::ranges::range Range> auto divisors(Range &&fac)
{
    using T = std::ranges::range_value_t<Range>::first_type;
    return divisors<T, Range>(std::forward<Range>(fac));
}

template <integral2 T, std::integral S = u32> std::vector<T> divisors(T num, const SPF<S> &spf = {})
{
    if (num == 1)
        return {1};
    return divisors(factor(std::move(num), spf));
}
} // namespace euler
