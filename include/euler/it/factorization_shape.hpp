#pragma once

#include "combinatorics.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates numbers with a given ordered factorization shape. For example, if shape is
/// `[a, b, c]`, this class will enumerate numbers of the form `p^a * q^b * r^c` with `p < q < r`
/// prime.
///
/// Usage: `nums_with_factorization_shape(primes, shape, limit)`.
///
/// Requirements:
/// * The range of primes is sorted ascending.
template <std::ranges::view Vp, std::ranges::view Vs, integral2 T>
class nums_with_ordered_factorization_shape : public it_base
{
  public:
    using value_type = T;

    nums_with_ordered_factorization_shape() = default;
    constexpr nums_with_ordered_factorization_shape(Vp primes, Vs shape, T limit)
        : _primes(std::move(primes)), _shape(std::move(shape)), _limit(limit)
    {
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _enumerate(std::ranges::begin(_primes), std::ranges::begin(_shape), T(1), f);
    }

    /// Reduces numbers with the given factorization shape. Uses std::transform_reduce.
    template <execution_policy Exec, typename U, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        auto its = std::ranges::begin(_shape);
        auto eTotal = std::accumulate(its, std::ranges::end(_shape), 0);
        assert(eTotal > 0);
        auto length = std::ranges::upper_bound(_primes, std::pow(_limit, 1.0 / eTotal)) - std::ranges::begin(_primes);
        return std::transform_reduce(std::forward<Exec>(exec), counting_iterator((int64_t)0), counting_iterator(length),
                                     std::move(init), std::forward<BinaryOp>(op), [&](auto i) {
                                         U total = init;
                                         _enumerate(std::ranges::begin(_primes) + i + 1, its + 1,
                                                    pow(T(_primes[i]), *its),
                                                    [&](auto &&x) { total = op(total, f(x)); });
                                         return total;
                                     });
    }

    /// Reduces numbers with the given factorization shape.
    template <typename U, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<U>)
    [[nodiscard]] U reduce(U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return it_base::reduce(std::move(init), std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
    }

  private:
    using Itp = std::ranges::iterator_t<const Vp>;
    using Its = std::ranges::iterator_t<const Vs>;
    using Tp = std::ranges::range_value_t<Vp>;
    using Ts = std::ranges::range_value_t<Vs>;

    Vp _primes;
    Vs _shape;
    T _limit;

    template <std::invocable<value_type> Fun> constexpr result_t _enumerate(Itp itp, Its its, T current, Fun f) const
    {
        if (its == _shape.end())
            return callbackResult(f, current);
        Ts e = *its;
        Ts eTotal = std::accumulate(its, _shape.end(), Ts(0));
        T bound = (T)std::pow(_limit / current, 1.0 / eTotal);
        for (; itp != std::ranges::end(_primes); ++itp)
        {
            Tp p = *itp;
            if (p > bound)
                break;
            if (!_enumerate(itp + 1, its + 1, current * pow(T(p), e), f))
                return it::result_break;
        }
        return it::result_continue;
    }
};

template <std::ranges::range Range1, std::ranges::range Range2, integral2 T>
nums_with_ordered_factorization_shape(Range1 &&, Range2 &&, T)
    -> nums_with_ordered_factorization_shape<std::views::all_t<Range1>, std::views::all_t<Range2>, T>;

/// Enumerates numbers with a given factorization shape up to a bound, using the given sorted list of primes. For
/// example, if shape is `[a, b, c]`, then this class will enumerate numbers of the form `p^a * q^b * r^c` with `p`,
/// `q`, `r` prime.
///
/// Usage: `nums_with_factorization_shape(primes, shape, limit)`.
///
/// Requirements:
/// * `primes` must be sorted ascending.
/// * `shape` must be sorted ascending.
template <std::ranges::view Vp, std::ranges::view Vs, integral2 T> class nums_with_factorization_shape : public it_base
{
  public:
    using value_type = T;

    nums_with_factorization_shape() = default;
    constexpr nums_with_factorization_shape(Vp primes, Vs shape, T limit)
        : _primes(std::move(primes)), _shape(std::move(shape)), _limit(limit)
    {
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return it::permutations(std::ranges::ref_view(_shape), std::ranges::size(_shape))([&](auto &&perm) {
            return it::nums_with_ordered_factorization_shape(
                std::ranges::ref_view(_primes), std::ranges::ref_view(perm), _limit)([&](auto &&n) { return f(n); });
        });
    }

    /// Reduces numbers with the given factorization shape. Uses std::transform_reduce.
    template <execution_policy Exec, typename U, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        U total = init;
        it::permutations(std::ranges::ref_view(_shape), std::ranges::size(_shape))([&](auto &&perm) {
            total = op(total, it::nums_with_ordered_factorization_shape(std::ranges::ref_view(_primes),
                                                                        std::ranges::ref_view(perm), _limit)
                                  .reduce(std::forward<Exec>(exec), std::move(init), std::forward<BinaryOp>(op),
                                          std::forward<UnaryOp>(f)));
        });
        return total;
    }

    /// Reduces numbers with the given factorization shape.
    template <typename U, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<U>)
    [[nodiscard]] U reduce(U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return it_base::reduce(std::move(init), std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
    }

  private:
    using Ts = std::ranges::range_value_t<Vs>;

    Vp _primes;
    Vs _shape;
    T _limit;
};

template <std::ranges::range Range1, std::ranges::range Range2, integral2 T>
nums_with_factorization_shape(Range1 &&, Range2 &&, T)
    -> nums_with_factorization_shape<std::views::all_t<Range1>, std::views::all_t<Range2>, T>;

/// Enumerates numbers with a given factorization shape up to a bound, using the given sorted list of primes. For
/// example, if shape is `[a, b, c]`, then this class will enumerate numbers of the form `p^a * q^b * r^c` with `p`,
/// `q`, `r` prime.
///
/// Usage: `nums_with_factorization_shape(primes, shape, limit)`.
///
/// Requirements:
/// * `primes` must be sorted ascending.
/// * `shape` must be sorted ascending.
///
/// This is an alternative implementation. It might be faster in some cases.
template <std::ranges::view Vp, std::ranges::view Vs, integral2 T> class nums_with_factorization_shape2 : public it_base
{
  public:
    using value_type = T;

    nums_with_factorization_shape2() = default;
    constexpr nums_with_factorization_shape2(Vp primes, Vs shape, T limit)
        : _primes(std::move(primes)), _shape(std::move(shape)), _limit(limit)
    {
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        std::vector<Tp> usedPrimes;
        return _enumerate(std::ranges::begin(_primes), std::ranges::begin(_shape), usedPrimes, T(1), f);
    }

    /// Reduces numbers with the given factorization shape. Uses std::transform_reduce.
    template <execution_policy Exec, typename U, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        if (std::ranges::empty(_shape))
            return U(1);
        auto its = std::ranges::begin(_shape);
        auto e = *its;
        auto bound = std::pow(_limit, 1.0 / e);
        auto stop = std::ranges::upper_bound(_primes, bound);
        return std::transform_reduce(std::forward<Exec>(exec), counting_iterator((int64_t)0),
                                     counting_iterator(stop - std::ranges::begin(_primes)), std::move(init),
                                     std::forward<BinaryOp>(op), [&](auto i) {
                                         auto p = _primes[i];
                                         T pe = pow(T(p), e);
                                         std::vector<Tp> usedPrimes{p};
                                         Itp itp2 = std::ranges::begin(_primes);
                                         // If next shape entry is same as previous, start from itp + 1.
                                         if (its + 1 != _shape.end() && *(its + 1) == e)
                                             itp2 += i + 1;
                                         U total{};
                                         _enumerate(itp2, its + 1, usedPrimes, pe, [&](auto &&x) { total += f(x); });
                                         return total;
                                     });
    }

    template <typename U, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<U>)
    [[nodiscard]] U reduce(U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return it_base::reduce(std::move(init), std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
    }

  private:
    using Itp = decltype(std::declval<Vp>().cbegin());
    using Its = decltype(std::declval<Vs>().cbegin());
    using Tp = std::ranges::range_value_t<Vp>;

    Vp _primes;
    Vs _shape;
    T _limit;

    template <std::invocable<value_type> Fun>
    constexpr result_t _enumerate(Itp itp, Its its, std::vector<Tp> &usedPrimes, T current, Fun f) const
    {
        if (its == _shape.end())
            return callbackResult(f, current);
        auto e = *its;
        auto bound = _limit / current;
        for (; itp != std::ranges::end(_primes); ++itp)
        {
            auto p = *itp;
            T pe = pow(T(p), e);
            if (pe > bound)
                break;
            if (std::ranges::contains(usedPrimes, p))
                continue;
            usedPrimes.push_back(p);
            Itp itp2 = std::ranges::begin(_primes);
            // If next shape entry is same as previous, start from itp + 1.
            if (its + 1 != _shape.end() && *(its + 1) == e)
                itp2 = itp + 1;
            if (!_enumerate(itp2, its + 1, usedPrimes, current * pe, f))
                return it::result_break;
            usedPrimes.pop_back();
        }
        return it::result_continue;
    }
};

template <std::ranges::range Range1, std::ranges::range Range2, integral2 T>
nums_with_factorization_shape2(Range1 &&, Range2 &&, T)
    -> nums_with_factorization_shape2<std::views::all_t<Range1>, std::views::all_t<Range2>, T>;
} // namespace it
} // namespace euler
