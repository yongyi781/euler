#pragma once

#include "../prime.hpp"
#include "factor.hpp"

namespace euler::it
{
/// Enumerates the factored divisors of a number.
template <integral2 T, std::ranges::view V> class divisors_factored_t : public it_base
{
    using It = std::ranges::iterator_t<const V>;
    V fac_;

  public:
    using prime_type = std::ranges::range_value_t<V>::first_type;
    using value_type = std::pair<T, const PF<prime_type> &>;

    divisors_factored_t() = default;
    constexpr divisors_factored_t(V fac) : fac_(std::move(fac)) {}

    template <typename Fun> result_t operator()(Fun f) const
    {
        PF<prime_type> fac_d;
        fac_d.data().reserve(std::ranges::size(fac_));
        auto const dfs = [&](this auto &&self, auto it, T cur) -> result_t {
            if (!callbackResult(f, value_type{cur, fac_d}))
                return result_break;
            for (; it != std::ranges::end(fac_); ++it)
            {
                auto &&[p, e] = *it;
                T n = cur * p;
                auto &entry = fac_d.emplace_back(p, 1);
                for (int i = 1; i <= e; ++i, ++entry.second, n *= p)
                    if (!self(std::ranges::next(it), n))
                        return result_break;
                fac_d.pop_back();
            }
            return result_continue;
        };
        return dfs(std::ranges::begin(fac_), T(1));
    }

    /// Returns the number of divisors, using a faster algorithm. Make sure to use a larger type `T`
    /// if you expect an answer larger than `2^64 - 1`.
    [[nodiscard]] constexpr T size() const
    {
        return it::wrap(std::ranges::ref_view(fac_)).map([](auto &&pe) -> T { return pe.second + 1; }).product();
    }
};

template <integral2 T, std::ranges::range Range>
it::divisors_factored_t<T, std::views::all_t<Range>> divisors_factored(Range &&fac)
{
    return {std::forward<Range>(fac)};
}

template <std::ranges::range Range>
it::divisors_factored_t<typename std::ranges::range_value_t<Range>::first_type, std::views::all_t<Range>>
divisors_factored(Range &&fac)
{
    return {std::forward<Range>(fac)};
}

template <integral2 T> it::divisors_factored_t<T, std::views::all_t<PF<T>>> divisors_factored(T num)
{
    return {PF{factor(std::move(num)).to()}};
}
} // namespace euler::it
