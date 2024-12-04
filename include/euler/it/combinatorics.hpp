#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates all possible combinations of `k` elements from a range.
template <std::ranges::view V> class combinations : public it_base
{
  public:
    using element_type = std::ranges::range_value_t<V>;
    using value_type = std::vector<element_type>;

    combinations() = default;
    constexpr combinations(V range, size_t k) : _range(std::move(range)), _k(k) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if (_k == 0)
            return callbackResult(f, value_type{});
        if (_k > std::ranges::size(_range))
            return result_continue;

        auto s = range(0UZ, _k - 1).map([&](size_t i) { return _range.begin() + i; }).to();
        while (true)
        {
            if (!it::callbackResult(f, s | std::views::transform([&](auto &&x) -> decltype(auto) { return *x; })))
                return it::result_break;
            auto i = _k - 1;
            while (++s[i] == _range.end() - _k + i + 1)
                if (--i == (size_t)-1)
                    return it::result_continue;
            for (auto j = i + 1; j < _k; ++j)
                s[j] = s[j - 1] + 1;
        }
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <typename Cont>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to() const
    {
        Cont result;
        (*this)([&](auto &&x) -> void {
            if constexpr (requires { result.emplace_back(x.begin(), x.end()); })
                result.emplace_back(x.begin(), x.end());
            else if constexpr (requires { result.push_back(value_type(x.begin(), x.end())); })
                result.push_back(value_type(x.begin(), x.end()));
            else if constexpr (requires { result.insert(value_type(x.begin(), x.end())); })
                result.insert(value_type(x.begin(), x.end()));
            else
                static_assert(false, "No acceptable insertion operation");
        });
        return result;
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <template <typename...> typename Cont = std::vector> [[nodiscard]] constexpr Cont<value_type> to() const
    {
        return to<Cont<value_type>>();
    }

    /// Collects the first n terms of this enumerable into a container.
    template <typename Cont>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to(size_t n) const
    {
        Cont result;
        result.reserve(n);
        size_t i = 0;
        (*this)([&](auto &&x) -> bool {
            if constexpr (requires { result.emplace_back(x.begin(), x.end()); })
                result.emplace_back(x.begin(), x.end());
            else if constexpr (requires { result.push_back(value_type(x.begin(), x.end())); })
                result.push_back(value_type(x.begin(), x.end()));
            else if constexpr (requires { result.insert(value_type(x.begin(), x.end())); })
                result.insert(value_type(x.begin(), x.end()));
            else
                static_assert(false, "No acceptable insertion operation");
            return ++i < n;
        });
        return result;
    }

    /// Collects the first n terms of this enumerable into a container.
    template <template <typename...> typename Cont = std::vector>
    [[nodiscard]] constexpr Cont<value_type> to(size_t n) const
    {
        return to<Cont<value_type>>(n);
    }

  private:
    V _range;
    size_t _k{};
};

template <std::ranges::range Range> combinations(Range &&, size_t) -> combinations<std::views::all_t<Range>>;

/// Enumerates all possible combinations of `k` elements from a range, with replacement.
template <std::ranges::view V> class combinations_with_replacement : public it_base
{
  public:
    using element_type = std::ranges::range_value_t<V>;
    using value_type = std::vector<element_type>;

    combinations_with_replacement() = default;
    constexpr combinations_with_replacement(V range, size_t k) : _range(std::move(range)), _k(k) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if (std::ranges::empty(_range))
            return result_continue;
        if (_k == 0)
            return callbackResult(f, value_type{});

        std::vector s(_k, std::ranges::begin(_range));
        while (true)
        {
            if (!callbackResult(f, s | std::views::transform([](auto &&x) -> decltype(auto) { return *x; })))
                return result_break;
            auto it = std::ranges::rbegin(s);
            while (++*it == std::ranges::end(_range))
            {
                if (++it == std::ranges::rend(s))
                    return result_continue;
                for (auto jt = it - 1; jt + 1 != std::ranges::rbegin(s); --jt)
                    *jt = *it + 1;
            }
        }
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <typename Cont>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to() const
    {
        Cont result;
        (*this)([&](auto &&x) -> void {
            if constexpr (requires { result.emplace_back(x.begin(), x.end()); })
                result.emplace_back(x.begin(), x.end());
            else if constexpr (requires { result.push_back(value_type(x.begin(), x.end())); })
                result.push_back(value_type(x.begin(), x.end()));
            else if constexpr (requires { result.insert(value_type(x.begin(), x.end())); })
                result.insert(value_type(x.begin(), x.end()));
            else
                static_assert(false, "No acceptable insertion operation");
        });
        return result;
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <template <typename...> typename Cont = std::vector> [[nodiscard]] constexpr Cont<value_type> to() const
    {
        return to<Cont<value_type>>();
    }

    /// Collects the first n terms of this enumerable into a container.
    template <typename Cont>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to(size_t n) const
    {
        Cont result;
        result.reserve(n);
        size_t i = 0;
        (*this)([&](auto &&x) -> bool {
            if constexpr (requires { result.emplace_back(x.begin(), x.end()); })
                result.emplace_back(x.begin(), x.end());
            else if constexpr (requires { result.push_back(value_type(x.begin(), x.end())); })
                result.push_back(value_type(x.begin(), x.end()));
            else if constexpr (requires { result.insert(value_type(x.begin(), x.end())); })
                result.insert(value_type(x.begin(), x.end()));
            else
                static_assert(false, "No acceptable insertion operation");
            return ++i < n;
        });
        return result;
    }

    /// Collects the first n terms of this enumerable into a container.
    template <template <typename...> typename Cont = std::vector>
    [[nodiscard]] constexpr Cont<value_type> to(size_t n) const
    {
        return to<Cont<value_type>>(n);
    }

  private:
    V _range;
    size_t _k{};
};

template <std::ranges::range Range>
combinations_with_replacement(Range &&, size_t) -> combinations_with_replacement<std::views::all_t<Range>>;

/// Enumerates all possible permutations of `k` elements from a range. The range must be sorted.
template <std::ranges::view V> class permutations : public it_base
{
  public:
    using value_type = combinations<V>::value_type;

    permutations() = default;
    constexpr permutations(V base, size_t k) : _combs(std::move(base), k) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        // Pass by value is intentional here.
        return _combs([&](auto &&comb) {
            std::vector v(std::ranges::begin(comb), std::ranges::end(comb));
            if (!callbackResult(f, v))
                return result_break;
            while (std::next_permutation(v.begin(), v.end()))
                if (!callbackResult(f, v))
                    return result_break;
            return result_continue;
        });
    }

  private:
    combinations<V> _combs;
};

template <std::ranges::range Range> permutations(Range &&, size_t) -> permutations<std::views::all_t<Range>>;

/// The `k`-fold self-Cartesian product of a range.
template <std::ranges::view V> class permutations_with_replacement : public it_base
{
  public:
    using value_type = combinations_with_replacement<V>::value_type;

    permutations_with_replacement() = default;
    constexpr permutations_with_replacement(V base, size_t k) : _combs(std::move(base), k) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _combs([&](auto &&comb) {
            std::vector v(std::ranges::begin(comb), std::ranges::end(comb));
            if (!callbackResult(f, v))
                return result_break;
            while (std::next_permutation(v.begin(), v.end()))
                if (!callbackResult(f, v))
                    return result_break;
            return result_continue;
        });
    }

  private:
    combinations_with_replacement<V> _combs;
};

template <std::ranges::range Range>
permutations_with_replacement(Range &&, size_t) -> permutations_with_replacement<std::views::all_t<Range>>;

/// The `k` fold self-Cartesian product of a range. This one is faster than non-sorted for smaller values of `k`.
template <std::ranges::view V> class power : public it_base
{
  public:
    using value_type = std::vector<std::ranges::range_value_t<V>>;

    power() = default;
    constexpr power(V base, size_t k) : _base(std::move(base)), _k(k) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if (_k == 0 || std::ranges::size(_base) == 1)
            return callbackResult(f, value_type(_k, *std::ranges::begin(_base)));
        if (_base.empty())
            return result_continue;
        std::vector s(_k, std::ranges::begin(_base));
        while (true)
        {
            if (!callbackResult(f, s | std::views::transform([](auto &&x) -> decltype(auto) { return *x; })))
                return result_break;
            auto i = s.rbegin();
            while (++*i == std::ranges::end(_base))
            {
                *i = std::ranges::begin(_base);
                if (++i == s.rend())
                    return result_continue;
            }
        }
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <typename Cont>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to() const
    {
        Cont result;
        (*this)([&](auto &&x) -> void {
            if constexpr (requires { result.emplace_back(x.begin(), x.end()); })
                result.emplace_back(x.begin(), x.end());
            else if constexpr (requires { result.push_back(value_type(x.begin(), x.end())); })
                result.push_back(value_type(x.begin(), x.end()));
            else if constexpr (requires { result.insert(value_type(x.begin(), x.end())); })
                result.insert(value_type(x.begin(), x.end()));
            else
                static_assert(false, "No acceptable insertion operation");
        });
        return result;
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <template <typename...> typename Cont = std::vector> [[nodiscard]] constexpr Cont<value_type> to() const
    {
        return to<Cont<value_type>>();
    }

    /// Collects the first n terms of this enumerable into a container.
    template <typename Cont>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to(size_t n) const
    {
        Cont result;
        result.reserve(n);
        size_t i = 0;
        (*this)([&](auto &&x) -> bool {
            if constexpr (requires { result.emplace_back(x.begin(), x.end()); })
                result.emplace_back(x.begin(), x.end());
            else if constexpr (requires { result.push_back(value_type(x.begin(), x.end())); })
                result.push_back(value_type(x.begin(), x.end()));
            else if constexpr (requires { result.insert(value_type(x.begin(), x.end())); })
                result.insert(value_type(x.begin(), x.end()));
            else
                static_assert(false, "No acceptable insertion operation");
            return ++i < n;
        });
        return result;
    }

    /// Collects the first n terms of this enumerable into a container.
    template <template <typename...> typename Cont = std::vector>
    [[nodiscard]] constexpr Cont<value_type> to(size_t n) const
    {
        return to<Cont<value_type>>(n);
    }

  private:
    V _base;
    size_t _k{};
};

template <std::ranges::range Range> power(Range &&, size_t) -> power<std::views::all_t<Range>>;

/// Cartesian product of a list of lists.
template <std::ranges::view V>
    requires std::ranges::range<std::ranges::range_value_t<V>>
class product : public it_base
{
  public:
    using value_type = std::vector<std::ranges::range_value_t<std::ranges::range_value_t<V>>>;

    product() = default;
    constexpr product(V items) : _items(std::move(items)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        value_type combination(_items.size());
        return _enumerate(combination, std::ranges::begin(_items), _items.end(), combination.begin(), f);
    }

  private:
    using It = std::ranges::iterator_t<const V>;
    using CIt = value_type::iterator;

    V _items;

    template <std::invocable<value_type> Fun>
    constexpr auto _enumerate(value_type &combination, It it, It end, CIt cit, Fun f) const
    {
        if (cit == combination.end())
            return callbackResult(f, combination);
        for (auto &&x : *it)
        {
            *cit = x;
            if (!_enumerate(combination, std::next(it), end, std::next(cit), f))
                return result_break;
        }
        return result_continue;
    }
};

/// Cartesian product of two enumerables
template <enumerable E1, enumerable E2> class product_it : public it_base
{
  public:
    using value_type = std::pair<typename E1::value_type, typename E2::value_type>;

    product_it() = default;
    constexpr product_it(E1 e1, E2 e2) : _e1(std::move(e1)), _e2(std::move(e2)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _e1([&](auto &&x) { return _e2([&](auto &&y) { return f(value_type{x, y}); }); });
    }

  private:
    E1 _e1;
    E2 _e2;
};

template <std::ranges::range Range> product(Range &&) -> product<std::views::all_t<Range>>;

/// Enumerates all subsets of a range. The range must have fewer than 64 elements.
template <std::ranges::view V> class powerset : public it_base
{
  public:
    using value_type = std::vector<std::ranges::range_value_t<V>>;

    powerset() = default;
    constexpr powerset(V range) : _range(std::move(range)) { assert(std::ranges::size(_range) < 64); }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        auto size = std::ranges::size(_range);
        value_type combination;
        combination.reserve(size);
        if (!callbackResult(f, combination))
            return result_break;
        for (uint64_t i = 1; i < 1ULL << size; ++i)
        {
            auto c = std::countr_zero(i);
            combination.erase(combination.end() - c, combination.end());
            combination.push_back(_range[size - c - 1]);
            if (!callbackResult(f, combination))
                return result_break;
        }
        return result_continue;
    }

  private:
    V _range;
};

template <std::ranges::range Range> powerset(Range &&) -> powerset<std::views::all_t<Range>>;

/// Enumerate nonnegative tuples with a given sum. Takes about 1 - 3.5 ns per tuple depending on length.
/// Usage: `it::tuples_with_sum(<length>, <total>)`.
template <integral2 T> class tuples_with_sum : public it_base
{
  public:
    using value_type = std::vector<T>;

    tuples_with_sum() = default;
    constexpr tuples_with_sum(size_t length, T total) : _length(length), _total(total) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        value_type t(_length, T(0));
        if (_length == 0)
        {
            if (_total == 0)
                return callbackResult(f, t);
            return result_break;
        }
        t.back() = _total;
        if (!callbackResult(f, t))
            return result_break;
        for (size_t index = t.size() - 2; index != (size_t)(-1); index--)
            if (!_enumerate(t, index, _total, f))
                return result_break;
        return result_continue;
    }

  private:
    size_t _length = 0;
    T _total = 0;

    template <std::invocable<value_type> Fun>
    constexpr result_t _enumerate(value_type &t, size_t index, T total, Fun f) const
    {
        for (T i = 1; i <= total; ++i)
        {
            t[index] = i;
            t.back() = total - i;
            if (!callbackResult(f, t))
                return result_break;
            if (total - i > 0)
                for (size_t index2 = t.size() - 2; index2 > index; index2--)
                    if (!_enumerate(t, index2, total - i, f))
                        return result_break;
            t[index] = 0;
        }
        return result_continue;
    }
};

/// Enumerates all partitions of `n`. Takes about 2.4 ns per partition generated.
template <std::integral T> class partitions : public it_base
{
  public:
    using value_type = std::vector<T>;

    partitions() = default;
    constexpr partitions(T n) : _n(n) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        value_type partition;
        partition.reserve(_n);
        if (!output(partition, _n, f))
            return result_break;
        return _enumerate(partition, _n, std::move(f));
    }

    /// Uses a faster algorithm to count of the number of partitions.
    template <integral2 Z = size_t> [[nodiscard]] constexpr Z size() const
    {
        std::vector<Z> p(_n + 1);
        p[0] = 1;
        for (T i = 1; i <= _n; ++i)
        {
            for (T k = 1;; ++k)
            {
                T j = (k * (3 * k - 1)) / 2;
                if (j > i)
                    break;
                if (k & 1)
                    p[i] += p[i - j];
                else
                    p[i] -= p[i - j];
                j = (k * (3 * k + 1)) / 2;
                if (j > i)
                    break;
                if (k & 1)
                    p[i] += p[i - j];
                else
                    p[i] -= p[i - j];
            }
        }
        return p[_n];
    }

  private:
    T _n;

    template <std::invocable<value_type> Fun> constexpr result_t output(value_type &partition, T n, Fun f) const
    {
        partition.push_back(n);
        auto res = callbackResult(f, partition);
        partition.pop_back();
        return res;
    }

    template <std::invocable<value_type> Fun> constexpr result_t _enumerate(value_type &partition, T n, Fun f) const
    {
        T lb = partition.empty() ? 1 : partition.back();
        for (T i = n / 2; i >= lb; --i)
        {
            partition.push_back(i);
            if (!output(partition, n - i, f))
                return result_break;
            if (i <= n / 3 && !_enumerate(partition, n - i, f))
                return result_break;
            partition.pop_back();
        }
        return result_continue;
    }
};
} // namespace it
} // namespace euler
