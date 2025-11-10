#pragma once

#include "../counting_iterator.hpp"
#include "../decls.hpp"
#include <boost/unordered/unordered_flat_set.hpp>
#include <cstdint>
#include <ranges>

/// @defgroup Stateless iterators
///
/// This file implements iterators via callbacks. Clang is very good at optimizing these. However, these iterators
/// cannot be paused, so some constructs are not possible, such as zipping two of these. (It is still possible zip one
/// of these with a classic iterator, though.)
namespace euler::it
{
enum result_t : uint8_t
{
    result_break,
    result_continue,
    result_continue_no_recurse
};

struct it_base;

template <typename T>
concept enumerable = std::derived_from<std::decay_t<T>, it_base>;

// Forward declarations.
template <enumerable E, std::invocable<typename E::value_type> Fn>
    requires(!std::is_void_v<std::invoke_result_t<Fn, typename E::value_type>>)
class map_t;

template <enumerable E> class take_t;

template <enumerable E> class cycle_t;

template <enumerable E, std::predicate<typename E::value_type> Pred> class take_while_t;

template <enumerable E> class drop_t;

template <enumerable E, std::predicate<typename E::value_type> Pred> class drop_while_t;

template <enumerable E, std::invocable<typename E::value_type> Fn>
    requires is_optional<std::invoke_result_t<Fn, typename E::value_type>>
class filter_map_t;

template <enumerable E, std::predicate<typename E::value_type> Pred> class filter_t;

template <enumerable E> class enumerate_t;

template <template <typename...> typename Set, enumerable E, std::invocable<typename E::value_type> Proj>
class unique_t;

/// Invoke a callable object, and returns result_continue if the callable returns void.
template <typename Callable, typename... Args> constexpr result_t callbackResult(Callable &&f, Args &&...args)
{
    if constexpr (std::is_void_v<std::invoke_result_t<Callable, Args...>>)
    {
        std::forward<Callable>(f)(std::forward<Args>(args)...);
        return result_continue;
    }
    else
    {
        return result_t(std::forward<Callable>(f)(std::forward<Args>(args)...));
    }
}

/// Base class for enumerables. `Tp` is the enumerable's value type.
struct it_base
{
    /// Gets the ith term of this enumerable. Returns default value if out of range.
    template <typename Self> [[nodiscard]] constexpr Self::value_type operator[](this const Self &self, size_t i)
    {
        typename Self::value_type result{};
        self([&](auto &&x) {
            if (i == 0)
            {
                result = std::forward<decltype(x)>(x);
                return result_break;
            }
            --i;
            return result_continue;
        });
        return result;
    }

    /// Gets the ith term of this enumerable, or nullopt if there isn't one.
    template <typename Self>
    [[nodiscard]] constexpr std::optional<typename Self::value_type> at(this const Self &self, size_t i)
    {
        std::optional<typename Self::value_type> result{};
        self([&](auto &&x) -> result_t {
            if (i == 0)
            {
                result = std::forward<decltype(x)>(x);
                return result_break;
            }
            --i;
            return result_continue;
        });
        return result;
    }

    /// Gets the number of elements in this enumerable. Note that this will run forever if the enumerable is infinite.
    template <typename Self> [[nodiscard]] constexpr size_t size(this const Self &self)
    {
        size_t result = 0;
        self([&](auto &&) -> void { ++result; });
        return result;
    }

    /// Gets the first term of this enumerable, or none if the enumerable is empty.
    template <typename Self>
    [[nodiscard]] constexpr std::optional<typename Self::value_type> first(this const Self &self)
    {
        std::optional<typename Self::value_type> result{};
        self([&](auto &&x) -> result_t {
            result = std::forward<decltype(x)>(x);
            return result_break;
        });
        return result;
    }

    /// Gets the last term of this enumerable, or none if the enumerable is empty. The enumerable must be finite.
    template <typename Self>
    [[nodiscard]] constexpr std::optional<typename Self::value_type> last(this const Self &self)
    {
        std::optional<typename Self::value_type> result{};
        self([&](auto &&x) -> void { result = std::forward<decltype(x)>(x); });
        return result;
    }

    /// Appends this enumerable into a container. The enumerable must be finite, or else this may run forever.
    template <typename Cont, typename Self>
        requires(!std::ranges::view<Cont>)
    constexpr void appendTo(this const Self &self, Cont &out)
    {
        self([&](auto &&x) -> void {
            if constexpr (requires { out.emplace_back(std::forward<decltype(x)>(x)); })
                out.emplace_back(std::forward<decltype(x)>(x));
            else if constexpr (requires { out.push_back(std::forward<decltype(x)>(x)); })
                out.push_back(std::forward<decltype(x)>(x));
            else
                out.insert(std::forward<decltype(x)>(x));
        });
    }

    /// Appends this enumerable into a container. The enumerable must be finite, or else this may run forever.
    template <typename Cont, typename Self>
        requires(!std::ranges::view<Cont>)
    constexpr void appendTo(this const Self &self, Cont &out, size_t n)
    {
        size_t i = 0;
        self([&](auto &&x) -> bool {
            if constexpr (requires { out.emplace_back(std::forward<decltype(x)>(x)); })
                out.emplace_back(std::forward<decltype(x)>(x));
            else if constexpr (requires { out.push_back(std::forward<decltype(x)>(x)); })
                out.push_back(std::forward<decltype(x)>(x));
            else
                out.insert(std::forward<decltype(x)>(x));
            return ++i < n;
        });
    }
    /// Collects this enumerable into a container. The enumerable must be finite.
    template <typename Cont, typename Self>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to(this const Self &self)
    {
        Cont result;
        self.appendTo(result);
        return result;
    }

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <template <typename...> typename Cont = std::vector, typename Self>
    [[nodiscard]] constexpr Cont<typename Self::value_type> to(this const Self &self)
    {
        return self.template to<Cont<typename Self::value_type>>();
    }

    /// Collects the first n terms of this enumerable into a container.
    template <typename Cont, typename Self>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to(this const Self &self, size_t n)
    {
        Cont result;
        result.reserve(n);
        self.appendTo(result, n);
        return result;
    }

    /// Collects the first n terms of this enumerable into a container.
    template <template <typename...> typename Cont = std::vector, typename Self>
    [[nodiscard]] constexpr Cont<typename Self::value_type> to(this const Self &self, size_t n)
    {
        return self.template to<Cont<typename Self::value_type>>(n);
    }

    /// Maps this enumerable by a function.
    template <typename Self, typename Fn> [[nodiscard]] constexpr map_t<Self, Fn> map(this Self self, Fn fn)
    {
        return {std::move(self), std::move(fn)};
    }

    /// Return elements from the enumerable until it is exhausted. Then repeat the sequence indefinitely.
    template <typename Self> [[nodiscard]] constexpr cycle_t<Self> cycle(this Self self) { return {std::move(self)}; }

    /// Takes the first n terms of this enumerable. This may move data from this enumerable.
    template <typename Self> [[nodiscard]] constexpr take_t<Self> take(this Self self, size_t n)
    {
        return {std::move(self), n};
    }

    /// Takes terms of this enumerable while the predicate is true.
    template <typename Self, typename Pred>
    [[nodiscard]] constexpr take_while_t<Self, Pred> takeWhile(this Self self, Pred pred)
    {
        return {std::move(self), std::move(pred)};
    }

    /// Drops the first n terms of this enumerable. This may move data from this enumerable.
    template <typename Self> [[nodiscard]] constexpr drop_t<Self> drop(this Self self, size_t n)
    {
        return drop_t{std::move(self), n};
    }

    /// Drops terms of this enumerable while the predicate is true.
    template <typename Self, typename Pred>
    [[nodiscard]] constexpr drop_while_t<Self, Pred> dropWhile(this Self self, Pred pred)
    {
        return {std::move(self), std::move(pred)};
    }

    /// Creates an iterator that both filters and maps.
    template <typename Self, typename Fn>
        requires is_optional<std::invoke_result_t<Fn, typename Self::value_type>>
    [[nodiscard]] constexpr filter_map_t<Self, Fn> filterMap(this Self self, Fn fn)
    {
        return {std::move(self), std::move(fn)};
    }

    /// Creates an iterator that filters.
    template <typename Self, typename Pred>
    [[nodiscard]] constexpr filter_t<Self, Pred> filter(this Self self, Pred pred)
    {
        return {std::move(self), std::move(pred)};
    }

    /// Numbers items from 0 onwards.
    template <typename Self> [[nodiscard]] constexpr enumerate_t<Self> enumerate(this Self self)
    {
        return {std::move(self)};
    }

    /// Outputs unique elements.
    template <template <typename...> typename Set = boost::unordered_flat_set, typename Self,
              typename Proj = std::identity>
    [[nodiscard]] constexpr unique_t<Set, Self, Proj> unique(this Self self, Proj proj = {})
    {
        return {std::move(self), std::move(proj)};
    }

    /// Reduces the enumerable into a single value.
    template <typename Self, typename T, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<T>)
    [[nodiscard]] constexpr T reduce(this const Self &self, T init, BinaryOp op, UnaryOp f = {})
    {
        T result = std::move(init);
        self([&](auto &&x) { result = op(std::move(result), f(std::forward<decltype(x)>(x))); });
        return result;
    }

    /// Calculates the sum of the enumerable.
    template <typename Self, execution_policy Exec, typename T = Self::value_type>
    [[nodiscard]] T sum(this const Self &self, Exec &&exec, T init = {})
    {
        return self.reduce(std::forward<Exec>(exec), std::move(init), std::plus{});
    }

    /// Calculates the sum of the enumerable.
    template <typename Self, typename T = Self::value_type>
        requires(!execution_policy<T>)
    [[nodiscard]] constexpr T sum(this const Self &self, T init = {})
    {
        return self.reduce(std::move(init), std::plus{});
    }

    /// Calculates the product of the enumerable. Uses std::transform_reduce.
    template <typename Self, execution_policy Exec, typename T = Self::value_type>
    [[nodiscard]] T product(this const Self &self, Exec &&exec, T init = T(1))
    {
        return self.reduce(std::forward<Exec>(exec), std::move(init), std::multiplies{});
    }

    /// Calculates the product of the enumerable.
    template <typename Self, typename T = Self::value_type>
        requires(!execution_policy<T>)
    [[nodiscard]] constexpr T product(this const Self &self, T init = T(1))
    {
        return self.reduce(std::move(init), std::multiplies{});
    }

    /// Calculates the minimum of the enumerable.
    template <typename Self, typename Proj = std::identity,
              std::strict_weak_order<std::invoke_result_t<Proj, typename Self::value_type>,
                                     std::invoke_result_t<Proj, typename Self::value_type>>
                  Comp = std::ranges::less>
    [[nodiscard]] constexpr Self::value_type min(this const Self &self, Comp comp = {}, Proj proj = {})
    {
        typename Self::value_type result{};
        bool init = false;
        self([&](auto &&x) -> void {
            if (!init || std::invoke(comp, std::invoke(proj, x), std::invoke(proj, result)))
            {
                init = true;
                result = x;
            }
        });
        return result;
    }

    /// Calculates the maximum of the enumerable.
    template <typename Self, typename Proj = std::identity,
              std::strict_weak_order<std::invoke_result_t<Proj, typename Self::value_type>,
                                     std::invoke_result_t<Proj, typename Self::value_type>>
                  Comp = std::ranges::less>
    [[nodiscard]] constexpr Self::value_type max(this const Self &self, Comp comp = {}, Proj proj = {})
    {
        typename Self::value_type result{};
        bool init = false;
        self([&](auto &&x) -> void {
            if (!init || std::invoke(comp, std::invoke(proj, result), std::invoke(proj, x)))
            {
                init = true;
                result = x;
            }
        });
        return result;
    }

    /// Returns true if the predicate is true for all items.
    template <typename Self, typename Pred = std::identity>
    [[nodiscard]] constexpr bool all(this const Self &self, Pred pred = {})
    {
        bool result = true;
        self([&](auto &&x) -> result_t {
            if (!std::invoke(pred, std::forward<decltype(x)>(x)))
            {
                result = false;
                return it::result_break;
            }
            return it::result_continue;
        });
        return result;
    }

    /// Returns true if the predicate is true for any item.
    template <typename Self, typename Pred = std::identity>
    [[nodiscard]] constexpr bool any(this const Self &self, Pred pred = {})
    {
        bool result = false;
        self([&](auto &&x) -> result_t {
            if (std::invoke(pred, std::forward<decltype(x)>(x)))
            {
                result = true;
                return it::result_break;
            }
            return it::result_continue;
        });
        return result;
    }
};

/// Maps this enumerable by a function.
template <enumerable E, std::invocable<typename E::value_type> Fn>
    requires(!std::is_void_v<std::invoke_result_t<Fn, typename E::value_type>>)
class map_t : public it_base
{
    E _base;
    Fn _fn;

  public:
    using value_type = std::remove_cvref_t<std::invoke_result_t<Fn, typename E::value_type>>;

    map_t() = default;
    constexpr map_t(E base, Fn fn) : _base(std::move(base)), _fn(std::move(fn)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        return _base([&](auto &&x) -> decltype(auto) { return f(std::invoke(_fn, std::forward<decltype(x)>(x))); });
    }

    /// Calculates the size by simply calling the size of the base enumerable.
    [[nodiscard]] constexpr size_t size() const { return std::ranges::size(_base); }

    /// Reduces the enumerable into a single value. Uses the base enumerable's reduce function.
    template <execution_policy Exec, typename T, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] T reduce(Exec &&exec, T init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return _base.reduce(std::forward<Exec>(exec), std::move(init), std::forward<BinaryOp>(op),
                            [&](auto &&x) { return f(std::invoke(_fn, std::forward<decltype(x)>(x))); });
    }

    /// Reduces the enumerable into a single value. Uses std::transform_reduce.
    template <typename T, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<T>)
    [[nodiscard]] constexpr T reduce(T init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return _base.reduce(std::move(init), std::forward<BinaryOp>(op),
                            [&](auto &&x) { return f(std::invoke(_fn, std::forward<decltype(x)>(x))); });
    }
};

/// Return elements from the enumerable until it is exhausted. Then repeat the sequence indefinitely.
template <enumerable E> class cycle_t : public it_base
{
    E _base;

  public:
    using value_type = E::value_type;

    cycle_t() = default;
    constexpr cycle_t(E base) : _base(std::move(base)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        while (true)
            if (!_base(f))
                break;
        return result_continue;
    }

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr size_t size() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto cycle() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto last() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto reduce() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto sum() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto product() const = delete;
};

template <std::ranges::range Range> cycle_t(Range &&) -> cycle_t<std::views::all_t<Range>>;

/// Takes the first `n` elements of this enumerable.
template <enumerable E> class take_t : public it_base
{
    E _base;
    size_t _n{};

  public:
    using value_type = E::value_type;

    take_t() = default;
    constexpr take_t(E base, size_t n) : _base(std::move(base)), _n(n) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        result_t result = result_continue;
        _base([&](auto &&x) {
            if (++i > _n)
                return result_break;
            return result = callbackResult(f, std::forward<decltype(x)>(x));
        });
        return result;
    }
};

/// Takes while the predicate is true.
template <enumerable E, std::predicate<typename E::value_type> Pred> class take_while_t : public it_base
{
    E _base;
    Pred _pred;

  public:
    using value_type = E::value_type;

    take_while_t() = default;
    constexpr take_while_t(E base, Pred pred) : _base(std::move(base)), _pred(std::move(pred)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        result_t result = result_continue;
        _base([&](auto &&x) {
            if (!_pred(x))
                return result_break;
            return result = callbackResult(f, std::forward<decltype(x)>(x));
        });
        return result;
    }
};

/// Drops the first `n` elements of this enumerable.
template <enumerable E> class drop_t : public it_base
{
    E _base;
    size_t _n{};

  public:
    using value_type = E::value_type;

    drop_t() = default;
    constexpr drop_t(E base, size_t n) : _base(std::move(base)), _n(n) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        return _base([&](auto &&x) {
            if (++i <= _n)
                return result_continue;
            return callbackResult(f, std::forward<decltype(x)>(x));
        });
    }
};

/// Drops while the predicate is true.
template <enumerable E, std::predicate<typename E::value_type> Pred> class drop_while_t : public it_base
{
    E _base;
    Pred _pred;

  public:
    using value_type = E::value_type;

    drop_while_t() = default;
    constexpr drop_while_t(E base, Pred pred) : _base(std::move(base)), _pred(std::move(pred)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        bool drop = true;
        return _base([&](auto &&x) {
            if (drop && _pred(x))
                return result_continue;
            drop = false;
            return callbackResult(f, std::forward<decltype(x)>(x));
        });
    }
};

/// An iterator that both filters and maps.
template <enumerable E, std::invocable<typename E::value_type> Fn>
    requires is_optional<std::invoke_result_t<Fn, typename E::value_type>>
class filter_map_t : public it_base
{
    E _base;
    Fn _fn;

  public:
    using value_type = std::invoke_result_t<Fn, typename E::value_type>::value_type;

    filter_map_t() = default;
    constexpr filter_map_t(E base, Fn fn) : _base(std::move(base)), _fn(std::move(fn)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        return _base([&](auto &&x) -> decltype(auto) {
            auto o = std::invoke(_fn, std::forward<decltype(x)>(x));
            return !o.has_value() || callbackResult(f, *o);
        });
    }
};

/// Take while the predicate is true.
template <enumerable E, std::predicate<typename E::value_type> Pred> class filter_t : public it_base
{
    E _base;
    Pred _pred;

  public:
    using value_type = E::value_type;

    filter_t() = default;
    constexpr filter_t(E base, Pred pred) : _base(std::move(base)), _pred(std::move(pred)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        return _base([&](auto &&x) {
            if (!_pred(x))
                return result_continue;
            return callbackResult(f, std::forward<decltype(x)>(x));
        });
    }
};

/// Numbers items from 0 onwards.
template <enumerable E> class enumerate_t : public it_base
{
    E _base;

  public:
    using value_type = std::pair<size_t, typename E::value_type>;

    enumerate_t() = default;
    constexpr enumerate_t(E e) : _base(std::move(e)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        return _base([&](auto &&x) { return callbackResult(f, value_type{i++, std::forward<decltype(x)>(x)}); });
    }

    /// This has the same size as the base enumerable.
    [[nodiscard]] constexpr size_t size() const { return std::ranges::size(_base); }
};

/// Outputs unique elements.
template <template <typename...> typename Set, enumerable E, std::invocable<typename E::value_type> Proj>
class unique_t : public it_base
{
    E _base;
    Proj _proj;

  public:
    using value_type = E::value_type;

    unique_t() = default;
    constexpr unique_t(E e, Proj proj) : _base(std::move(e)), _proj(std::move(proj)) {}

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        using T = std::remove_cvref_t<std::invoke_result_t<Proj, value_type>>;
        Set<T> visited;
        return _base([&](auto &&x) -> result_t {
            auto key = std::invoke(_proj, std::forward<decltype(x)>(x));
            if (!visited.contains(key))
            {
                visited.insert(key);
                return callbackResult(f, std::forward<decltype(x)>(x));
            }
            return it::result_continue;
        });
    }
};

/// Wraps a `std::ranges::view` as an enumerable.
template <std::ranges::view V> class wrap : public it_base
{
    V _base;

  public:
    using value_type = std::ranges::range_value_t<V>;

    wrap() = default;
    constexpr wrap(V base) : _base(std::move(base)) {}

    /// Invokes f with the specified execution policy. Note that in this variant, the enumeration
    /// cannot be broken out of.
    template <execution_policy Exec, typename Fun> result_t operator()(Exec &&exec, Fun f) const
    {
        std::for_each(std::forward<Exec>(exec), std::ranges::begin(_base), std::ranges::end(_base), std::move(f));
        return result_continue;
    }

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        for (auto &&x : _base)
            if (!callbackResult(f, std::forward<decltype(x)>(x)))
                return result_break;
        return result_continue;
    }

    /// Reduces the enumerable into a single value. Uses std::transform_reduce.
    template <execution_policy Exec, typename T, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] T reduce(Exec &&exec, T init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return std::transform_reduce(std::forward<Exec>(exec), std::ranges::begin(_base), std::ranges::end(_base),
                                     std::move(init), std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
    }

    /// Reduces the enumerable into a single value. Uses std::transform_reduce.
    template <typename T, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<T>)
    [[nodiscard]] constexpr T reduce(T init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        return std::transform_reduce(std::ranges::begin(_base), std::ranges::end(_base), std::move(init),
                                     std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
    }
};

template <std::ranges::range Range> wrap(Range &&) -> wrap<std::views::all_t<Range>>;

/// Iterates over the interval `[begin, begin + step, ..., end]`.
template <integral2 T> class range : public it_base
{
    T _begin;
    T _end;
    T _step;

  public:
    using value_type = T;

    range() = default;
    constexpr range(T begin, T end, T step = T(1))
        : _begin(std::move(begin)), _end(std::move(end)), _step(std::move(step))
    {
    }

    /// Invokes `f` with the specified execution policy. Note that in this variant, the enumeration cannot be broken out
    /// of.
    template <execution_policy Exec, typename Fun> result_t operator()(Exec &&exec, Fun f) const
    {
        if (!empty())
        {
            if (_step == 1)
            {
                std::for_each(std::forward<Exec>(exec), counting_iterator{_begin}, counting_iterator{T(_end + 1)},
                              std::move(f));
            }
            else
            {
                T n = (_end - _begin) / _step;
                std::for_each(std::forward<Exec>(exec), counting_iterator{T(0)}, counting_iterator{T(n + 1)},
                              [&](auto &&i) { f(i * _step + _begin); });
            }
        }
        return result_continue;
    }

    template <typename Fun> constexpr result_t operator()(Fun f) const
    {
        if (_step > 0)
        {
            for (T i = _begin; i <= _end; i += _step)
                if (!callbackResult(f, i))
                    return result_break;
        }
        else
        {
            for (T i = _begin; i >= _end; i += _step)
                if (!callbackResult(f, i))
                    return result_break;
        }
        return result_continue;
    }

    /// Returns whether this range is empty.
    [[nodiscard]] constexpr bool empty() const { return _step == 0 || (_step > 0 ? _end < _begin : _begin < _end); }

    [[nodiscard]] constexpr size_t size() const { return empty() ? 0 : (_end - _begin) / _step + 1; }

    /// Reduces the enumerable into a single value. Uses `std::transform_reduce`.
    template <execution_policy Exec, typename U, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] U reduce(Exec &&exec, U init, BinaryOp op, UnaryOp f = {}) const
    {
        if (empty())
            return init;
        if (_step == 1)
            return std::transform_reduce(std::forward<Exec>(exec), counting_iterator{_begin},
                                         counting_iterator{T(_end + 1)}, std::move(init), std::move(op), std::move(f));
        T n = (_end - _begin) / _step;
        return std::transform_reduce(std::forward<Exec>(exec), counting_iterator{T(0)}, counting_iterator{T(n + 1)},
                                     std::move(init), std::move(op), [&](auto &&i) { return f(i * _step + _begin); });
    }

    /// Reduces the enumerable into a single value. Uses `std::transform_reduce`.
    template <typename U, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<U>)
    [[nodiscard]] constexpr U reduce(U init, BinaryOp op, UnaryOp f = {}) const
    {
        if (empty())
            return init;
        if (_step == 1)
            return std::transform_reduce(counting_iterator{_begin}, counting_iterator{T(_end + 1)}, std::move(init),
                                         std::move(op), std::move(f));
        return it_base::reduce(std::move(init), std::move(op), std::move(f));
    }

    /// Returns the first element for which `pred` returns `true`.
    template <typename Pred> [[nodiscard]] std::optional<T> findIf(Pred pred) const
    {
        if (empty())
            return std::nullopt;
        if (_step == 1)
        {
            auto const res = std::find_if(counting_iterator{_begin}, counting_iterator{T(_end + 1)}, std::move(pred));
            if (res.value == _end + 1)
                return std::nullopt;
            return res.value;
        }
        T n = (_end - _begin) / _step;
        auto const res = std::find_if(counting_iterator{T(0)}, counting_iterator{T(n + 1)},
                                      [&](auto &&i) { return pred(i * _step + _begin); });
        if (res.value == n + 1)
            return std::nullopt;
        return res.value * _step + _begin;
    }

    /// Returns, in parallel, the first element for which `pred` returns `true`.
    template <execution_policy Exec, typename Pred> [[nodiscard]] std::optional<T> findIf(Exec &&exec, Pred pred) const
    {
        if (empty())
            return std::nullopt;
        if (_step == 1)
        {
            auto const res = std::find_if(std::forward<Exec>(exec), counting_iterator{_begin},
                                          counting_iterator{T(_end + 1)}, std::move(pred));
            if (res.value == _end + 1)
                return std::nullopt;
            return res.value;
        }
        T n = (_end - _begin) / _step;
        auto const res = std::find_if(std::forward<Exec>(exec), counting_iterator{T(0)}, counting_iterator{T(n + 1)},
                                      [&](auto &&i) { return pred(i * _step + _begin); });
        if (res.value == n + 1)
            return std::nullopt;
        return res.value * _step + _begin;
    }
};

template <integral2 TBegin, integral2 TEnd> range(TBegin, TEnd) -> range<std::common_type_t<TBegin, TEnd>>;

template <integral2 TBegin, integral2 TEnd, integral2 TStep>
range(TBegin, TEnd, TStep) -> range<std::common_type_t<TBegin, TEnd, TStep>>;
// Some helper functions

/// Prints an enumerable to a stream.
template <typename CharT, typename Traits, it::enumerable E>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const E &e)
{
    bool isFirst = true;
    e([&](auto &&x) {
        if (isFirst)
            isFirst = false;
        else
            o << ", ";
        o << x;
    });
    return o;
}
} // namespace euler::it
