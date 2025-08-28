#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates the Farey sequence of rationals with a bounded denominator between 0 and 1, inclusive.
template <integral2 T> class farey_sequence : public it::it_base
{
    T _n;

  public:
    using value_type = std::pair<T, T>;

    constexpr farey_sequence(T n) : _n(std::move(n)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        T a = 0;
        T b = 1;
        T c = 1;
        T d = _n;
        while (a <= _n)
        {
            if (!callbackResult(f, std::pair{a, b}))
                return result_break;
            T k = (_n + b) / d;
            std::tie(a, b, c, d) = std::tuple{c, d, k * c - a, k * d - b};
        }
        return it::result_continue;
    }
};
} // namespace it
} // namespace euler
