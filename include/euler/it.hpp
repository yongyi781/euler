#pragma once

#include "decls.hpp"
#include <boost/unordered/unordered_flat_set_fwd.hpp>
#include <execution>
#include <fstream>
#include <ranges>
#include <utility>

inline namespace euler
{
/// @defgroup Stateless iterators
///
/// This file implements iterators via callbacks. Clang is very good at optimizing these. However, these iterators
/// cannot be paused, so some constructs are not possible, such as zipping two of these. (It is still possible zip one
/// of these with a classic iterator, though.)
namespace it
{
enum result_t : uint8_t
{
    result_break,
    result_continue,
    result_continue_no_recurse
};

/// Predicate that always returns true.
constexpr auto pred_true = [](auto &&) { return true; };

// Forward declarations.
template <enumerable E, std::invocable<typename E::value_type> Fn> class map_t;

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
template <typename Callable, typename... Args>
    requires std::invocable<Callable, Args...>
constexpr result_t callbackResult(Callable &&f, Args &&...args) noexcept(std::is_nothrow_invocable_v<Callable, Args...>)
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

    /// Collects this enumerable into a container. The enumerable must be finite.
    template <typename Cont, typename Self>
        requires(!std::ranges::view<Cont>)
    [[nodiscard]] constexpr Cont to(this const Self &self)
    {
        Cont result;
        self([&](auto &&x) -> void {
            if constexpr (requires { result.emplace_back(std::forward<decltype(x)>(x)); })
                result.emplace_back(std::forward<decltype(x)>(x));
            else if constexpr (requires { result.push_back(std::forward<decltype(x)>(x)); })
                result.push_back(std::forward<decltype(x)>(x));
            else
                result.insert(std::forward<decltype(x)>(x));
        });
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
        size_t i = 0;
        self([&](auto &&x) -> bool {
            if constexpr (requires { result.emplace_back(std::forward<decltype(x)>(x)); })
                result.emplace_back(std::forward<decltype(x)>(x));
            else if constexpr (requires { result.push_back(std::forward<decltype(x)>(x)); })
                result.push_back(std::forward<decltype(x)>(x));
            else
                result.insert(std::forward<decltype(x)>(x));
            return ++i < n;
        });
        return result;
    }

    /// Collects the first n terms of this enumerable into a container.
    template <template <typename...> typename Cont = std::vector, typename Self>
    [[nodiscard]] constexpr Cont<typename Self::value_type> to(this const Self &self, size_t n)
    {
        return self.template to<Cont<typename Self::value_type>>(n);
    }

    /// Maps this enumerable by a function.
    template <typename Self, std::invocable<typename Self::value_type> Fn>
    [[nodiscard]] constexpr map_t<Self, Fn> map(this Self self, Fn fn)
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
    template <typename Self, std::predicate<typename Self::value_type> Pred>
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
    template <typename Self, std::predicate<typename Self::value_type> Pred>
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
    template <typename Self, std::predicate<typename Self::value_type> Pred>
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
              std::invocable<typename Self::value_type> Proj = std::identity>
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
    [[nodiscard]] constexpr T sum(this const Self &self, Exec &&exec, T init = {})
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
    [[nodiscard]] constexpr T product(this const Self &self, Exec &&exec, T init = T(1))
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
    template <typename Self, std::invocable<typename Self::value_type> Proj = std::identity,
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
    template <typename Self, std::invocable<typename Self::value_type> Proj = std::identity,
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
    template <typename Self, std::predicate<typename Self::value_type> Pred = std::identity>
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
    template <typename Self, std::predicate<typename Self::value_type> Pred = std::identity>
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

/// Wraps a `std::ranges::view` as an enumerable.
template <std::ranges::view V> class wrap : public it_base
{
  public:
    using value_type = std::ranges::range_value_t<V>;

    wrap() = default;
    constexpr wrap(V base) : _base(std::move(base)) {}

    /// Invokes f with the specified execution policy. Note that in this variant, the enumeration
    /// cannot be broken out of.
    template <execution_policy Exec, std::invocable<value_type> Fun>
    constexpr result_t operator()(Exec &&exec, Fun f) const
    {
        std::for_each(std::forward<Exec>(exec), std::ranges::begin(_base), std::ranges::end(_base), std::move(f));
        return result_continue;
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        for (auto &&x : _base)
            if (!callbackResult(f, std::forward<decltype(x)>(x)))
                return result_break;
        return result_continue;
    }

    /// Reduces the enumerable into a single value. Uses std::transform_reduce.
    template <execution_policy Exec, typename T, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] constexpr T reduce(Exec &&exec, T init, BinaryOp &&op, UnaryOp &&f = {}) const
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

  private:
    V _base;
};

template <std::ranges::range Range> wrap(Range &&) -> wrap<std::views::all_t<Range>>;

/// Iterates over the interval `[begin, begin + step, ..., end]`.
template <integral2 T> class range : public it_base
{
  public:
    using value_type = T;

    range() = default;
    constexpr range(T begin, T end, T step = T(1))
        : _begin(std::move(begin)), _end(std::move(end)), _step(std::move(step))
    {
    }

    /// Invokes `f` with the specified execution policy. Note that in this variant, the enumeration cannot be broken out
    /// of.
    template <execution_policy Exec, std::invocable<value_type> Fun>
    constexpr result_t operator()(Exec &&exec, Fun f) const
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

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        for (T i = _begin; _step > 0 ? i <= _end : i >= _end; i += _step)
            if (!callbackResult(f, i))
                return result_break;
        return result_continue;
    }

    /// Returns whether this range is empty.
    [[nodiscard]] constexpr bool empty() const { return _step == 0 || (_step > 0 ? _end < _begin : _begin < _end); }

    [[nodiscard]] constexpr size_t size() const { return empty() ? 0 : (_end - _begin) / _step + 1; }

    /// Reduces the enumerable into a single value. Uses `std::transform_reduce`.
    template <execution_policy Exec, typename U, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] constexpr U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        if (empty())
            return init;
        if (_step == 1)
            return std::transform_reduce(std::forward<Exec>(exec), counting_iterator{_begin},
                                         counting_iterator{T(_end + 1)}, std::move(init), std::forward<BinaryOp>(op),
                                         std::forward<UnaryOp>(f));
        T n = (_end - _begin) / _step;
        return std::transform_reduce(std::forward<Exec>(exec), counting_iterator{T(0)}, counting_iterator{T(n + 1)},
                                     std::move(init), std::forward<BinaryOp>(op),
                                     [&](auto &&i) { return f(i * _step + _begin); });
    }

    /// Reduces the enumerable into a single value. Uses `std::transform_reduce`.
    template <typename U, typename BinaryOp, typename UnaryOp = std::identity>
        requires(!execution_policy<U>)
    [[nodiscard]] constexpr U reduce(U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        if (empty())
            return init;
        if (_step == 1)
            return std::transform_reduce(counting_iterator{_begin}, counting_iterator{T(_end + 1)}, std::move(init),
                                         std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
        return it_base::reduce(std::move(init), std::forward<BinaryOp>(op), std::forward<UnaryOp>(f));
    }

  private:
    T _begin;
    T _end;
    T _step;
};

template <integral2 TBegin, integral2 TEnd> range(TBegin, TEnd) -> range<std::common_type_t<TBegin, TEnd>>;

template <integral2 TBegin, integral2 TEnd, integral2 TStep>
range(TBegin, TEnd, TStep) -> range<std::common_type_t<TBegin, TEnd, TStep>>;

/// Maps this enumerable by a function.
template <enumerable E, std::invocable<typename E::value_type> Fn> class map_t : public it_base
{
  public:
    using value_type = std::remove_cvref_t<std::invoke_result_t<Fn, typename E::value_type>>;

    map_t() = default;
    constexpr map_t(E base, Fn fn) : _base(std::move(base)), _fn(std::move(fn)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _base([&](auto &&x) -> decltype(auto) { return f(std::invoke(_fn, std::forward<decltype(x)>(x))); });
    }

    /// Calculates the size by simply calling the size of the base enumerable.
    [[nodiscard]] constexpr size_t size() const { return std::ranges::size(_base); }

    /// Reduces the enumerable into a single value. Uses the base enumerable's reduce function.
    template <execution_policy Exec, typename T, typename BinaryOp, typename UnaryOp = std::identity>
    [[nodiscard]] constexpr T reduce(Exec &&exec, T init, BinaryOp &&op, UnaryOp &&f = {}) const
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

  private:
    E _base;
    Fn _fn;
};

/// Return elements from the enumerable until it is exhausted. Then repeat the sequence indefinitely.
template <enumerable E> class cycle_t : public it_base
{
  public:
    using value_type = E::value_type;

    cycle_t() = default;
    constexpr cycle_t(E base) : _base(std::move(base)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
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

  private:
    E _base;
};

template <std::ranges::range Range> cycle_t(Range &&) -> cycle_t<std::views::all_t<Range>>;

/// Takes the first `n` elements of this enumerable.
template <enumerable E> class take_t : public it_base
{
  public:
    using value_type = E::value_type;

    take_t() = default;
    constexpr take_t(E base, size_t n) : _base(std::move(base)), _n(n) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
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

  private:
    E _base;
    size_t _n{};
};

/// Takes while the predicate is true.
template <enumerable E, std::predicate<typename E::value_type> Pred> class take_while_t : public it_base
{
  public:
    using value_type = E::value_type;

    take_while_t() = default;
    constexpr take_while_t(E base, Pred pred) : _base(std::move(base)), _pred(std::move(pred)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        result_t result = result_continue;
        _base([&](auto &&x) {
            if (!_pred(x))
                return result_break;
            return result = callbackResult(f, std::forward<decltype(x)>(x));
        });
        return result;
    }

  private:
    E _base;
    Pred _pred;
};

/// Drops the first `n` elements of this enumerable.
template <enumerable E> class drop_t : public it_base
{
  public:
    using value_type = E::value_type;

    drop_t() = default;
    constexpr drop_t(E base, size_t n) : _base(std::move(base)), _n(n) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        return _base([&](auto &&x) {
            if (++i <= _n)
                return result_continue;
            return callbackResult(f, std::forward<decltype(x)>(x));
        });
    }

  private:
    E _base;
    size_t _n{};
};

/// Drops while the predicate is true.
template <enumerable E, std::predicate<typename E::value_type> Pred> class drop_while_t : public it_base
{
  public:
    using value_type = E::value_type;

    drop_while_t() = default;
    constexpr drop_while_t(E base, Pred pred) : _base(std::move(base)), _pred(std::move(pred)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        bool drop = true;
        return _base([&](auto &&x) {
            if (drop && _pred(x))
                return result_continue;
            drop = false;
            return callbackResult(f, std::forward<decltype(x)>(x));
        });
    }

  private:
    E _base;
    Pred _pred;
};

/// An iterator that both filters and maps.
template <enumerable E, std::invocable<typename E::value_type> Fn>
    requires is_optional<std::invoke_result_t<Fn, typename E::value_type>>
class filter_map_t : public it_base
{
  public:
    using value_type = std::invoke_result_t<Fn, typename E::value_type>::value_type;

    filter_map_t() = default;
    constexpr filter_map_t(E base, Fn fn) : _base(std::move(base)), _fn(std::move(fn)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _base([&](auto &&x) -> decltype(auto) {
            auto o = std::invoke(_fn, std::forward<decltype(x)>(x));
            return !o.has_value() || callbackResult(f, *o);
        });
    }

  private:
    E _base;
    Fn _fn;
};

/// Take while the predicate is true.
template <enumerable E, std::predicate<typename E::value_type> Pred> class filter_t : public it_base
{
  public:
    using value_type = E::value_type;

    filter_t() = default;
    constexpr filter_t(E base, Pred pred) : _base(std::move(base)), _pred(std::move(pred)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _base([&](auto &&x) {
            if (!_pred(x))
                return result_continue;
            return callbackResult(f, std::forward<decltype(x)>(x));
        });
    }

  private:
    E _base;
    Pred _pred;
};

/// Numbers items from 0 onwards.
template <enumerable E> class enumerate_t : public it_base
{
  public:
    using value_type = std::pair<size_t, typename E::value_type>;

    enumerate_t() = default;
    constexpr enumerate_t(E e) : _base(std::move(e)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        return _base([&](auto &&x) { return callbackResult(f, value_type{i++, std::forward<decltype(x)>(x)}); });
    }

    /// This has the same size as the base enumerable.
    [[nodiscard]] constexpr size_t size() const { return std::ranges::size(_base); }

  private:
    E _base;
};

/// Outputs unique elements.
template <template <typename...> typename Set, enumerable E, std::invocable<typename E::value_type> Proj>
class unique_t : public it_base
{
  public:
    using value_type = E::value_type;

    unique_t() = default;
    constexpr unique_t(E e, Proj proj) : _base(std::move(e)), _proj(std::move(proj)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
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

  private:
    E _base;
    Proj _proj;
};

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

/// Enumerates a virtual tree. This version enumerates iteratively, rather than recursively.
/// Usage: `it::tree(root, fun(x, f, <call f on each child of x>), fun(x, <whether x has children>))`.
template <typename T, typename Fun, std::predicate<T> Pred = decltype(pred_true)> class tree : public it_base
{
  public:
    using value_type = T;

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    /// @param hasChildren A function that returns whether a node has children. Optional, but really speeds up the
    ///                    runtime if provided.
    constexpr tree(T root, Fun childrenFun, Pred hasChildren)
        : _root(std::move(root)), _childrenFun(std::move(childrenFun)), _hasChildren(std::move(hasChildren))
    {
    }

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    constexpr tree(T root, Fun childrenFun) : tree(std::move(root), std::move(childrenFun), pred_true) {}

    template <std::invocable<value_type> Callback> constexpr result_t operator()(Callback f) const
    {
        if (!it::callbackResult(f, _root))
            return result_break;
        std::vector<T> s{std::move(_root)};
        while (!s.empty())
        {
            T current = std::move(s.back());
            s.pop_back();
            if (!it::callbackResult(_childrenFun, std::move(current), [&](T child) -> result_t {
                    if (!it::callbackResult(f, child))
                        return result_break;
                    if (std::invoke(_hasChildren, child))
                        s.emplace_back(std::move(child));
                    return result_continue;
                }))
                return result_break;
        }
        return result_continue;
    }

  private:
    T _root;
    Fun _childrenFun;
    Pred _hasChildren;
};

/// Enumerates a virtual tree. This version enumerates recursively and in preorder traversal order.
/// Usage: `it::tree_preorder(root, fun(x, f, <call f on each child of x>), fun(x, <whether x has children>))`.
template <typename T, typename Fun, std::predicate<T> Pred = decltype(pred_true)> class tree_preorder : public it_base
{
  public:
    using value_type = T;

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    /// @param hasChildren A function that returns whether a node has children. Optional, but really speeds up the
    ///                    runtime if provided.
    constexpr tree_preorder(T root, Fun childrenFun, Pred hasChildren)
        : _root(std::move(root)), _childrenFun(std::move(childrenFun)), _hasChildren(std::move(hasChildren))
    {
    }

    /// Constructs an instance of this class.
    /// @param root The root of the tree.
    /// @param childrenFun A function that enumerates the children of a node. If you want inorder traversal, this
    ///                    should enumerate children backwards.
    constexpr tree_preorder(T root, Fun childrenFun) : tree_preorder(std::move(root), std::move(childrenFun), pred_true)
    {
    }

    template <std::invocable<value_type> Callback> constexpr result_t operator()(Callback f) const
    {
        if (!callbackResult(f, _root))
            return result_break;
        return dfs(_root, std::move(f));
    }

  private:
    T _root;
    Fun _childrenFun;
    Pred _hasChildren;

    template <std::invocable<T> Callback> result_t dfs(T node, Callback f) const
    {
        return callbackResult(_childrenFun, std::move(node), [&](T child) -> result_t {
            if (!callbackResult(f, child))
                return result_break;
            if (!std::invoke(_hasChildren, child))
                return result_continue;
            return dfs(std::move(child), f);
        });
    }
};

/// Enumerates smooth numbers up to a given limit.
template <std::ranges::view V, integral2 T> class smooth_numbers : public it_base
{
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

  private:
    using It = std::ranges::iterator_t<const V>;

    V _primes;
    T _limit;
};

/// Enumerates smooth numbers up to a given limit, along with their factorizations.
/// Usage: `smooth_numbers_factored(<primes>, <limit>)`.
template <std::ranges::view V, integral2 T> class smooth_numbers_factored : public it_base
{
  public:
    using value_type = std::pair<T, Factorization<T>>;
    using value_reference_type = std::pair<T, const Factorization<T> &>;

    smooth_numbers_factored() = default;
    constexpr smooth_numbers_factored(V primes, T limit) : _primes(std::move(primes)), _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        Factorization<T> fac{};
        fac.reserve(4);
        return it::tree_preorder(
            std::tuple<T, Factorization<T> &, It>{T(1), fac, std::ranges::begin(_primes)},
            [&](auto &&x, auto f) {
                auto &&[n, fac, i] = x;
                bool exists = !fac.empty();
                for (auto j = i; j != std::ranges::end(_primes) && mulLeq(n, *j, _limit); ++j)
                {
                    if (exists)
                        ++fac.back().second;
                    else
                        fac.emplace_back(*j, 1);
                    if (!f({n * *j, fac, j}))
                        return result_break;
                    if (exists)
                    {
                        --fac.back().second;
                        exists = false;
                    }
                    else
                        fac.pop_back();
                }
                return result_continue;
            },
            [&](auto &&x) { return mulLeq(get<0>(x), *get<2>(x), _limit); })(
            [&](auto &&x) { return f(value_reference_type{get<0>(x), get<1>(x)}); });
    }

  private:
    using It = std::ranges::iterator_t<const V>;

    V _primes;
    T _limit;
};

template <std::ranges::range Range, integral2 T>
smooth_numbers_factored(Range &&, T) -> smooth_numbers_factored<std::views::all_t<Range>, T>;

template <std::ranges::range Range, integral2 T>
smooth_numbers(Range &&, T) -> smooth_numbers<std::views::all_t<Range>, T>;

/// Enumerates squarefree smooth numbers up to a given limit.
template <std::ranges::view V, integral2 T> class squarefree_smooth_numbers : public it_base
{
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

  private:
    using It = std::ranges::iterator_t<const V>;

    V _primes;
    T _limit;
};

template <std::ranges::range Range, integral2 T>
squarefree_smooth_numbers(Range &&, T) -> squarefree_smooth_numbers<std::views::all_t<Range>, T>;

template <integral2 T = int64_t> class fibonacci : public it_base
{
  public:
    using value_type = T;

    constexpr fibonacci(T a = 0, T b = 1) : _a(a), _b(b) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        T x = _a;
        T y = _b;

        while (true)
        {
            if (!callbackResult(f, x))
                return result_break;
            std::tie(x, y) = std::pair{y, T(x + y)};
        }
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

  private:
    T _a, _b;
};

/// Enumerates primitive Pythaogrean triples. Takes about 7 ns per triple. Fast fast!
template <integral2 T = int64_t> class primitive_pythagorean_triples : public it_base
{
  public:
    using value_type = std::array<T, 3>;

    primitive_pythagorean_triples() = default;
    constexpr primitive_pythagorean_triples(T limit) : _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _enumerate(T(3), T(4), T(5), f);
    }

    [[nodiscard]] constexpr T limit() const { return _limit; }

  private:
    T _limit;

    template <std::invocable<value_type> Fun> constexpr result_t _enumerate(T a, T b, T c, Fun f) const
    {
        if (c > _limit)
            return result_continue;
        result_t result = callbackResult(f, value_type{a, b, c});
        if (result == result_break)
            return result_break;
        if (result == result_continue_no_recurse)
            return result_continue;
        std::array arr{std::array{2 * a + b - c, -2 * a + 2 * b + 2 * c, -2 * a + b + 3 * c},
                       std::array{2 * a + b + c, 2 * a - 2 * b + 2 * c, 2 * a - b + 3 * c},
                       std::array{2 * a - b + c, 2 * a + 2 * b + 2 * c, 2 * a + b + 3 * c}};
        for (auto &&[x, y, z] : arr)
        {
            result = _enumerate(x, y, z, f);
            if (result == result_break)
                return result_break;
            if (result == result_continue_no_recurse)
                return result_continue;
        }
        return result_continue;
    }
};

/// Enumerates primitive Pythaogrean triples. Takes about 27 ns per triple.
template <integral2 T = int64_t> class primitive_pythagorean_triples2 : public it_base
{
  public:
    using value_type = std::array<T, 3>;

    primitive_pythagorean_triples2() = default;
    constexpr primitive_pythagorean_triples2(T limit) : _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        using std::gcd;
        using euler::gcd;

        T hv = isqrt(_limit / 2);
        for (T v = 1; v <= hv; ++v)
        {
            for (T u = v + 1; u * u + v * v <= _limit; u += 2)
            {
                if (gcd(u, v) != 1)
                    continue;
                if (!callbackResult(f, value_type{u * u - v * v, 2 * u * v, u * u + v * v}))
                    return result_break;
            }
        }
        return result_continue;
    }

    constexpr T limit() const { return _limit; }

  private:
    T _limit;
};

/// Enumerates Pythaogrean triples. Takes about 0.5 ns per triple if using the faster version
/// of the primitive triple generator. Fast fast!
template <typename T = int64_t, typename P = primitive_pythagorean_triples<T>>
class pythagorean_triples : public it_base
{
  public:
    using value_type = P::value_type;

    pythagorean_triples() = default;
    constexpr pythagorean_triples(T limit) : base(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return base([&](auto &&x) {
            auto &&[a, b, c] = std::forward<decltype(x)>(x);
            T hk = base.limit() / c;
            for (T k = 1; k <= hk; ++k)
                if (!callbackResult(f, value_type{k * a, k * b, k * c}))
                    return result_break;
            return result_continue;
        });
    }

  private:
    P base;
};

/// Enumerates the digits of a number in a specified base, which must be ≥ 2.
template <integral2 T = int64_t, integral2 TBase = int> class digits : public it_base
{
  public:
    using value_type = TBase;

    digits() = default;
    constexpr digits(T n, TBase base = 10) : _n(n), _base(base) { assert(base >= 2); }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if constexpr (boost::integer_traits<T>::digits > 128)
            if (_base == 10)
            {
                for (char c : to_string(_n) | std::views::reverse)
                    if (!callbackResult(f, TBase(c - '0')))
                        return result_break;
                return result_continue;
            }
        T n = _n;
        while (n > 0)
        {
            if (!callbackResult(f, TBase(n % _base)))
                return result_break;
            n /= _base;
        }
        return result_continue;
    }

  private:
    T _n;
    TBase _base;
};

/// Unfolds an infinite recursive sequence starting from `a₀` using the function `fn`.
///
/// Usage: `unfold(a₀, fn)` where `fn: T → T` or `fn: T × ℕ → T`.
///
/// For example, `a₁` will be either `fn(a₀)` or `fn(a₀, 1)`.
template <typename T, typename Fn>
    requires(std::invocable<Fn, T> || std::invocable<Fn, T, size_t>)
class unfold : public it_base
{
  public:
    using value_type = T;

    unfold() = default;
    constexpr unfold(T a0, Fn fn) : _a0(std::move(a0)), _fn(std::move(fn)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        T state = _a0;
        while (true)
        {
            if (!callbackResult(f, state))
                return result_break;
            if constexpr (std::invocable<Fn, T, size_t>)
            {
                if constexpr (std::same_as<std::invoke_result_t<Fn, T, std::size_t>, std::optional<T>>)
                {
                    auto result = _fn(state, ++i);
                    if (!result)
                        return result_break;
                    state = *result;
                }
                else
                    state = _fn(state, ++i);
            }
            else if constexpr (std::same_as<std::invoke_result_t<Fn, T>, std::optional<T>>)
            {
                auto result = _fn(state);
                if (!result)
                    return result_break;
                state = *result;
            }
            else
            {
                state = _fn(state);
            }
        }
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

  private:
    T _a0;
    Fn _fn;
};

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
    [[nodiscard]] constexpr U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        auto its = std::ranges::begin(_shape);
        auto eTotal = std::accumulate(its, std::ranges::end(_shape), 0);
        assert(eTotal > 0);
        auto length = std::ranges::upper_bound(_primes, std::pow(_limit, 1.0 / eTotal)) - std::ranges::begin(_primes);
        return std::transform_reduce(std::forward<Exec>(exec), counting_iterator(0LL), counting_iterator(length),
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
    [[nodiscard]] constexpr U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
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
    [[nodiscard]] constexpr U reduce(Exec &&exec, U init, BinaryOp &&op, UnaryOp &&f = {}) const
    {
        if (std::ranges::empty(_shape))
            return U(1);
        auto its = std::ranges::begin(_shape);
        auto e = *its;
        auto bound = std::pow(_limit, 1.0 / e);
        auto stop = std::ranges::upper_bound(_primes, bound);
        return std::transform_reduce(std::forward<Exec>(exec), counting_iterator(0LL),
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

/// Enumerates all product partitions of `n`.
template <integral2 T> class product_partitions : public it_base
{
  public:
    using value_type = std::vector<T>;

    product_partitions() = default;
    constexpr product_partitions(T n) : _n(std::move(n)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        value_type partition;
        return _enumerate(partition, _n, f);
    }

  private:
    T _n;

    template <std::invocable<value_type> Fun> constexpr result_t _enumerate(value_type &partition, T n, Fun f) const
    {
        partition.push_back(n);
        if (!callbackResult(f, partition))
            return result_break;
        partition.pop_back();

        T lb = partition.empty() ? T(2) : partition.back();
        for (T i = isqrt(n); i >= lb; i--)
        {
            if (n % i == 0)
            {
                partition.push_back(i);
                if (!_enumerate(partition, n / i, f))
                    return result_break;
                partition.pop_back();
            }
        }
        return result_continue;
    }
};

/// Enumerates the Farey sequence of rationals with a bounded denominator between 0 and 1, inclusive.
template <integral2 T> class farey_sequence : public it::it_base
{
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

  private:
    T _n;
};

/// Splits a string at a character. TODO: make this work for wstrings too.
class split : public it_base
{
  public:
    using value_type = std::string;

    split(std::string s, char d) : _s(std::move(s)), _d(d) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        std::istringstream iss(_s);
        for (std::string token; std::getline(iss, token, _d);)
            if (!it::callbackResult(f, std::move(token)))
                return it::result_break;
        return it::result_continue;
    }

  private:
    std::string _s;
    char _d;
};

/// Gets lines from a file.
class lines : public it_base
{
  public:
    using value_type = std::string;

    lines(std::string fileName) : _fileName(std::move(fileName)) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        std::ifstream fin(_fileName);
        for (std::string line; std::getline(fin, line);)
            if (!it::callbackResult(f, std::move(line)))
                return it::result_break;
        return it::result_continue;
    }

  private:
    std::string _fileName;
};
} // namespace it

// Some helper functions

/// Prints an enumerable to a stream.
template <typename CharT, typename Traits, enumerable E>
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

template <integral2 TBase = int>
constexpr std::vector<TBase> digits(const integral2 auto &num, TBase base = 10,
                                    std::endian endian = std::endian::little)
{
    std::vector<TBase> result;
    result.reserve(8); // Easy optimization
    it::digits(num, base)([&](TBase d) { result.push_back(d); });
    if (endian == std::endian::big)
        reverse(result.begin(), result.end());
    return result;
}

template <integral2 T, std::ranges::range Range>
constexpr T fromDigits(Range &&digits, int base = 10, std::endian endian = std::endian::little)
{
    T result = 0;
    T p = 1;
    if (endian == std::endian::big)
        for (auto it = digits.rbegin(); it != digits.rend(); ++it)
        {
            result += p * *it;
            p *= base;
        }
    else
        for (auto it = digits.begin(); it != digits.end(); ++it)
        {
            result += p * *it;
            p *= base;
        }
    return result;
}

/// Shortcut for `it::digits(num, base).size()`.
template <integral2 T, integral2 TBase = int> constexpr size_t countDigits(const T &num, TBase base = 10)
{
    return it::digits(num, base).size();
}

/// Sums the digits of a number.
template <integral2 T, integral2 TBase = int> constexpr TBase sumDigits(const T &num, TBase base = 10)
{
    return it::digits(num, base).sum();
}

/// Returns the number with the reversed digits of n.
template <integral2 T, integral2 TBase = int> constexpr T reverseDigits(const T &n, TBase base = 10)
{
    return it::digits(n, base).reduce(T(0), [&](auto &&a, auto &&b) -> T { return base * a + b; });
}

/// Returns whether a number is a palindrome.
template <integral2 T, integral2 TBase = int> constexpr bool isPalindrome(const T &n, TBase base = 10)
{
    return n == reverseDigits(n, base);
}

/// Calculates the digit signature of a given number.
template <int Base = 10> constexpr std::array<int, Base> getDigitSignature(const integral2 auto &n)
{
    std::array<int, Base> result{};
    it::digits(n, Base)([&](int d) { ++result[d]; });
    return result;
}
} // namespace euler
