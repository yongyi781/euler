#pragma once

#include "../PF.hpp"
#include "tree.hpp"

namespace euler
{
namespace it
{
/// Enumerates smooth numbers up to a given limit.
template <std::ranges::view V, integral2 T> class smooth_numbers : public it_base
{
    using It = std::ranges::iterator_t<const V>;

    V _primes;
    T _limit;

  public:
    using value_type = T;

    smooth_numbers() = default;
    constexpr smooth_numbers(V primes, T limit) : _primes(std::move(primes)), _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return tree(
            std::pair{T(1), std::ranges::begin(_primes)},
            [&](auto &&x, auto f) {
                auto &&[n, i] = x;
                for (auto j = i; j != std::ranges::end(_primes) && mulLeq(n, *j, _limit); ++j)
                    if (!f(std::pair{n * *j, j}))
                        return it::result_break;
                return it::result_continue;
            },
            [&](auto &&x) { return mulLeq(x.first, *x.second, _limit); })([&](auto &&x) { return f(x.first); });
    }
};

template <std::ranges::range Range, integral2 T>
smooth_numbers(Range &&, T) -> smooth_numbers<std::views::all_t<Range>, T>;

/// Enumerates smooth numbers up to a given limit, along with their factorizations.
/// Usage: `smooth_numbers_factored(<primes>, <limit>)`.
template <std::ranges::view V, integral2 T> class smooth_numbers_factored : public it_base
{
    using It = std::ranges::iterator_t<const V>;

    V _primes;
    T _limit;

  public:
    using value_type = std::pair<T, PF<T>>;
    using value_reference_type = std::pair<T, const PF<T> &>;

    smooth_numbers_factored() = default;
    constexpr smooth_numbers_factored(V primes, T limit) : _primes(std::move(primes)), _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        PF<T> fac{};
        fac.data().reserve(4);
        return it::tree_preorder(
            std::tuple<T, PF<T> &, It>{T(1), fac, std::ranges::begin(_primes)},
            [&](auto &&x, auto f) {
                auto &&[n, fac, i] = x;
                bool exists = !fac.empty();
                for (auto j = i; j != std::ranges::end(_primes) && mulLeq(n, *j, _limit); ++j)
                {
                    if (exists)
                        ++fac.back().second;
                    else
                        fac.data().emplace_back(*j, 1);
                    if (!f({n * *j, fac, j}))
                        return result_break;
                    if (exists)
                    {
                        --fac.back().second;
                        exists = false;
                    }
                    else
                        fac.data().pop_back();
                }
                return result_continue;
            },
            [&](auto &&x) { return mulLeq(get<0>(x), *get<2>(x), _limit); })(
            [&](auto &&x) { return f(value_reference_type{get<0>(x), get<1>(x)}); });
    }
};

template <std::ranges::range Range, integral2 T>
smooth_numbers_factored(Range &&, T) -> smooth_numbers_factored<std::views::all_t<Range>, T>;

/// Enumerates squarefree smooth numbers up to a given limit.
template <std::ranges::view V, integral2 T> class squarefree_smooth_numbers : public it_base
{
    using It = std::ranges::iterator_t<const V>;

    V _primes;
    T _limit;

  public:
    using value_type = T;

    squarefree_smooth_numbers() = default;
    constexpr squarefree_smooth_numbers(V primes, T limit) : _primes(std::move(primes)), _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return tree(
            std::pair{T(1), std::ranges::begin(_primes)},
            [&](auto &&x, auto f) {
                auto &&[n, i] = x;
                for (auto j = i; j != std::ranges::end(_primes) && mulLeq(n, *j, _limit); ++j)
                    if (!f(std::pair{n * *j, j + 1}))
                        return it::result_break;
                return it::result_continue;
            },
            [&](auto &&x) { return x.second != std::ranges::end(_primes) && mulLeq(x.first, *x.second, _limit); })(
            [&](auto &&x) { return f(x.first); });
    }
};

template <std::ranges::range Range, integral2 T>
squarefree_smooth_numbers(Range &&, T) -> squarefree_smooth_numbers<std::views::all_t<Range>, T>;
} // namespace it
} // namespace euler
