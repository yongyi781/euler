#pragma once

#include "decls.hpp"
#include <map>
#include <random>
#include <set>

inline namespace euler
{
/// Returns the range [begin, end], possibly transformed by a function if given.
template <integral2 T, integral2 U, std::invocable<std::common_type_t<T, U>> Fun = std::identity>
constexpr auto range(T begin, U end, Fun f = {})
{
    using Tp = std::common_type_t<T, U>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, Tp>>;

    Tp b = begin;
    Tp e = end;
    if (e < b)
        return std::vector<FT>{};
    std::vector<FT> result((size_t)(e - b + 1));
    auto it = result.begin();
    for (Tp i = b; i <= e; ++i)
        *it++ = f(i);
    return result;
}

/// Returns the range [begin, end], possibly transformed by a function if given.
template <integral2 T, integral2 U, integral2 V, std::invocable<std::common_type_t<T, U, V>> Fun = std::identity>
constexpr auto range(T begin, U end, V step, Fun f = {})
{
    using Tp = std::common_type_t<T, U, V>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, Tp>>;

    Tp b = begin;
    Tp e = end;
    Tp s = step;

    if (s == 0 || (s > 0 && e < b) || (s < 0 && e > b))
        return std::vector<FT>{};
    std::vector<FT> result((size_t)((e - b) / s + 1));
    auto it = result.begin();
    if (s > 0)
        for (Tp i = b; i <= e; i += s)
            *it++ = f(i);
    else
        for (Tp i = b; i >= e; i += s)
            *it++ = f(i);
    return result;
}

/// Returns an array of values in the range [Begin, End], possibly transformed by a function if
/// given. This gives a performance advantage since it performs no allocations.
template <int64_t Begin, int64_t End, std::invocable<int64_t> Fun = std::identity>
    requires(End >= Begin)
constexpr auto arange(Fun f = {})
{
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, int64_t>>;
    std::array<FT, End - Begin + 1> result;
    auto it = result.begin();
    for (int64_t i = Begin; i <= End; ++i)
        *it++ = f(i);
    return result;
}

/// Returns an array of values in the range [Begin, End], possibly transformed by a function if
/// given.
template <int64_t Begin, int64_t End, int64_t Step, std::invocable<int64_t> Fun = std::identity>
    requires(Step != 0 && (Step > 0 ? End >= Begin : End <= Begin))
constexpr auto arange(Fun f = {})
{
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, int64_t>>;
    std::array<FT, (size_t)((End - Begin) / Step + 1)> result;
    auto it = result.begin();
    if constexpr (Step > 0)
        for (int64_t i = Begin; i <= End; i += Step)
            *it++ = f(i);
    else
        for (int64_t i = Begin; i >= End; i += Step)
            *it++ = f(i);
    return result;
}

/// Unfolds a recursive sequence into a vector. The parameter `f` can take either the previous
/// element of the sequence or the pair (previous element, index) as parameter.
template <typename T, typename Fun>
    requires std::invocable<Fun, T, size_t> || std::invocable<Fun, T>
constexpr std::vector<T> unfold(size_t size, T seed, Fun f)
{
    std::vector<T> result(size);
    result[0] = seed;
    for (size_t i = 1; i < size; ++i)
        if constexpr (std::invocable<Fun, T, size_t>)
            result[i] = f(result[i - 1], i);
        else
            result[i] = f(result[i - 1]);
    return result;
}

/// Transforms a range.
template <execution_policy Exec, std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> Fun>
constexpr auto mapv(Exec &&exec, Range &&r, Fun f)
{
    using T = std::ranges::range_value_t<Range>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
    std::vector<FT> result(r.size());
    std::transform(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), std::ranges::begin(result),
                   std::move(f));
    return result;
}

/// Transforms a range.
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> Fun>
constexpr auto mapv(Range &&r, Fun f)
{
    using T = std::ranges::range_value_t<Range>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
    std::vector<FT> result(r.size());
    std::transform(std::ranges::begin(r), std::ranges::end(r), std::ranges::begin(result), std::move(f));
    return result;
}

/// Reduces a range.
template <execution_policy Exec, std::ranges::range Range, typename T,
          std::invocable<std::ranges::range_value_t<Range>> Fun = std::identity,
          std::invocable<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>,
                         std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>
              BinaryOp>
constexpr T reducev(Exec &&exec, Range &&r, T init, BinaryOp op, Fun f = {})
{
    return std::transform_reduce(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), std::move(init),
                                 std::move(op), std::move(f));
}

/// Reduces a range.
template <std::ranges::range Range, typename T, std::invocable<std::ranges::range_value_t<Range>> Fun = std::identity,
          std::invocable<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>,
                         std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>
              BinaryOp>
constexpr T reducev(Range &&r, T init, BinaryOp op, Fun f = {})
{
    return std::transform_reduce(std::ranges::begin(r), std::ranges::end(r), std::move(init), std::move(op),
                                 std::move(f));
}

/// Reduces a range.
template <execution_policy Exec, integral2 T, integral2 U, typename Tp,
          std::invocable<std::common_type_t<T, U>> Fun = std::identity,
          std::invocable<std::invoke_result_t<Fun, std::common_type_t<T, U>>,
                         std::invoke_result_t<Fun, std::common_type_t<T, U>>>
              BinaryOp>
constexpr Tp reduceRange(Exec &&exec, T begin, U end, Tp init, BinaryOp op, Fun f = {})
{
    using V = std::common_type_t<T, U>;
    if (V(end) < V(begin))
        return init;
    return std::transform_reduce(std::forward<Exec>(exec), counting_iterator(V(begin)), counting_iterator(V(end + 1)),
                                 std::move(init), std::move(op), std::move(f));
}

/// Reduces a range.
template <integral2 T, integral2 U, typename Tp, std::invocable<std::common_type_t<T, U>> Fun = std::identity,
          std::invocable<std::invoke_result_t<Fun, std::common_type_t<T, U>>,
                         std::invoke_result_t<Fun, std::common_type_t<T, U>>>
              BinaryOp>
constexpr Tp reduceRange(T begin, U end, Tp init, BinaryOp op, Fun f = {})
{
    using V = std::common_type_t<T, U>;
    if (V(end) < V(begin))
        return init;
    return std::transform_reduce(counting_iterator(V(begin)), counting_iterator(V(end + 1)), std::move(init),
                                 std::move(op), std::move(f));
}

/// Sums a function over a range.
template <execution_policy Exec, std::ranges::range Range,
          std::invocable<std::ranges::range_value_t<Range>> Fun = std::identity>
auto sum(Exec &&exec, Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Exec>(exec), std::forward<Range>(r), T{}, std::plus{}, std::move(f));
}

/// Sums a function over a range.
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> Fun = std::identity>
constexpr auto sum(Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Range>(r), T{}, std::plus{}, std::move(f));
}

/// Sums a function over a range of numbers.
template <execution_policy Exec, integral2 T, integral2 U, std::invocable<std::common_type_t<T, U>> Fun = std::identity>
auto sum(Exec &&exec, T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end), Tp{}, std::plus{}, std::move(f));
}

/// Sums a function over a range of numbers.
template <integral2 T, integral2 U, std::invocable<std::common_type_t<T, U>> Fun = std::identity>
constexpr auto sum(T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::move(begin), std::move(end), Tp{}, std::plus{}, std::move(f));
}

/// Multiplies a function over a range.
template <execution_policy Exec, std::ranges::range Range,
          std::invocable<std::ranges::range_value_t<Range>> Fun = std::identity>
auto product(Exec &&exec, Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Exec>(exec), std::forward<Range>(r), T{1}, std::multiplies{}, std::move(f));
}

/// Multiplies a function over a range.
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> Fun = std::identity>
constexpr auto product(Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Range>(r), T{1}, std::multiplies{}, std::move(f));
}

/// Multiplies a function over a range of numbers.
template <execution_policy Exec, integral2 T, integral2 U, std::invocable<std::common_type_t<T, U>> Fun = std::identity>
auto product(Exec &&exec, T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end), Tp{1}, std::multiplies{},
                       std::move(f));
}

/// Multiplies a function over a range of numbers.
template <integral2 T, integral2 U, std::invocable<std::common_type_t<T, U>> Fun = std::identity>
constexpr auto product(T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::move(begin), std::move(end), Tp{1}, std::multiplies{}, std::move(f));
}

/// Returns the maximum value in a sequence according to a specified key selector function.
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> Key = std::identity>
constexpr auto maxBy(Range &&range, std::ranges::range_value_t<Range> init, Key key = {})
{
    return reducev(range, std::move(init), [&](auto &&a, auto &&b) {
        return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
    });
}

/// Returns the maximum value in a sequence according to a specified key selector function.
template <execution_policy Exec, std::ranges::range Range,
          std::invocable<std::ranges::range_value_t<Range>> Key = std::identity>
constexpr auto maxBy(Exec &&exec, Range &&range, std::ranges::range_value_t<Range> init, Key key = {})
{
    return reducev(std::forward<Exec>(exec), std::forward<Range>(range), std::move(init), [&](auto &&a, auto &&b) {
        return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
    });
}

/// Returns the maximum value in a numeric range according to a specified key selector function.
template <integral2 T, integral2 U, typename Tp, std::invocable<std::common_type_t<T, U>> Key = std::identity>
constexpr auto maxRangeBy(T begin, U end, Tp init, Key key = {})
{
    return reduceRange(std::move(begin), std::move(end), std::move(init), [&](auto &&a, auto &&b) {
        return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
    });
}

/// Returns the maximum value in a numeric range according to a specified key selector function.
template <execution_policy Exec, integral2 TBegin, integral2 TEnd, typename Tp,
          std::invocable<std::common_type_t<TBegin, TEnd>> Key = std::identity>
constexpr auto maxRangeBy(Exec &&exec, TBegin begin, TEnd end, Tp init, Key key = {})
{
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end), std::move(init),
                       [&](auto &&a, auto &&b) {
                           return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
                       });
}

/// Generates a random vector of a given size.
template <typename T, typename Gen = std::mt19937_64>
    requires std::integral<T> || std::floating_point<T>
std::vector<T> randomVector(size_t size, T low, T high, Gen &&g = {})
{
    std::vector<T> result(size);
    if constexpr (std::integral<T>)
    {
        std::uniform_int_distribution dist(low, high);
        std::generate(result.begin(), result.end(), [&]() { return dist(g); });
    }
    else
    {
        std::uniform_real_distribution dist(low, high);
        std::generate(result.begin(), result.end(), [&]() { return dist(g); });
    }
    return result;
}

/// Returns the least integer `n` such that `f(n) ≥ target`.
template <integral2 T, integral2 U, std::invocable<U> Fun>
constexpr U bisectionLowerBound(Fun f, const T &target, U low, U high)
{
    while (high - low > 1)
    {
        U mid = low + (high - low) / 2;
        if (f(mid) < target)
            low = mid;
        else
            high = mid;
    }
    return high;
}

/// Returns the least integer `n` such that `f(n) > target`.
template <integral2 T, integral2 U, std::invocable<U> Fun>
constexpr U bisectionUpperBound(Fun f, const T &target, U low, U high)
{
    while (high - low > 1)
    {
        U mid = low + (high - low) / 2;
        if (f(mid) > target)
            high = mid;
        else
            low = mid;
    }
    return high;
}

/// Finds a root of the function `f`.
template <typename T, std::invocable<T> Fun> constexpr T findRoot(Fun f, T low, T high, size_t maxSteps = 1000)
{
    auto fl = f(low);
    auto fh = f(high);
    if (fl == 0)
        return low;
    if (fh == 0)
        return high;
    if ((fl > 0) == (fh > 0))
        throw std::domain_error("Function has same sign at endpoints");
    if (fl > 0)
        std::swap(low, high);
    while (maxSteps-- > 0)
    {
        T mid = (low + high) / 2;
        if (low == mid || high == mid || f(mid) == 0)
            return mid;
        if (f(mid) < 0)
            low = mid;
        else
            high = mid;
    }
    return high;
}

template <std::ranges::range Range> constexpr void flatten(Range &&v, std::ranges::range auto &out)
{
    if constexpr (std::ranges::range<std::ranges::range_value_t<Range>>)
        for (const auto &e : std::forward<Range>(v))
            flatten(e, out);
    else
        out.insert(out.end(), std::ranges::begin(v), std::ranges::end(v));
}

/// Returns the (recursive) total size of a range of ranges.
template <std::ranges::sized_range Range> constexpr size_t totalSize(Range &&v)
{
    using T = std::ranges::range_value_t<Range>;
    if constexpr (std::ranges::sized_range<T>)
        return sum(std::forward<Range>(v), [](auto &&x) { return totalSize(x); });
    else
        return std::ranges::size(std::forward<Range>(v));
}

/// Flattens a range of ranges into a single vector.
template <typename T, std::ranges::range Range> constexpr std::vector<T> flatten(Range &&v)
{
    std::vector<T> result;
    result.reserve(totalSize(v));
    for (const auto &e : std::forward<Range>(v))
        flatten(e, result);
    return result;
}

/// Same as `collections.Counter` in Python.
template <template <typename...> typename Map = std::map, std::ranges::range Range> auto counter(Range &&r)
{
    using T = std::ranges::range_value_t<Range>;
    Map<T, size_t> result;
    for (auto &&i : std::forward<Range>(r))
        ++result[i];
    return result;
}

// Period finding

struct period_result
{
    size_t period;
    size_t preperiod;

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const period_result &x)
    {
        std::basic_ostringstream<CharT, Traits> ss;
        ss.flags(o.flags());
        ss.imbue(o.getloc());
        ss.precision(o.precision());
        ss << "{ period: " << x.period << ", preperiod: " << x.preperiod << " }";
        return o << std::move(ss).str();
    }
};

/// Finds the period of a recurrence given by `a(0) = a0, a(n) = f(a(n-1)).`.
template <typename T, std::invocable<T> Fun> period_result findPeriod(Fun f, T a0, double ratio = 2.0)
{
    size_t periodBound = 1;
    period_result res{.period = 1};
    T tortoise = a0;
    T hare = f(a0);
    while (tortoise != hare)
    {
        if (res.period == periodBound)
        {
            tortoise = hare;
            periodBound = std::max(periodBound + 1, (size_t)((double)periodBound * ratio));
            res.period = 0;
        }
        hare = f(hare);
        ++res.period;
    }
    tortoise = a0;
    hare = a0;
    for (size_t i = 0; i < res.period; ++i)
        hare = f(hare);
    while (tortoise != hare)
    {
        tortoise = f(tortoise);
        hare = f(hare);
        ++res.preperiod;
    }
    return res;
}

/// Sums a known periodic function by taking advantage of the period.
template <typename Fun, integral2 Z> auto sumPeriodic(Fun fn, Z preperiod, Z period, Z start, Z stop)
{
    using T = decltype(fn(Z(0)));
    T total = sum(start, std::min(stop, preperiod - 1), ([&](auto &&i) -> decltype(auto) { return fn(i); }));
    T mid = 0;
    start = std::max(preperiod, start);
    Z const l = (stop - start + 1) % period + start - 1;
    for (Z i = start; i <= std::min(stop, start + period - 1); ++i)
    {
        auto res = fn(i);
        mid += res;
        if (i <= l)
            total += res;
    }
    return total + (stop - start + 1) / period * mid;
}

/// Appends a range to a vector.
template <typename T, std::ranges::range Range> constexpr std::vector<T> &operator+=(std::vector<T> &a, Range &&b)
{
    a.insert(a.end(), std::ranges::begin(b), std::ranges::end(b));
    return a;
}

/// Appends a range to a set.
template <typename T, std::ranges::range Range> constexpr std::set<T> &operator+=(std::set<T> &a, Range &&b)
{
    a.insert(std::ranges::begin(b), std::ranges::end(b));
    return a;
}

/// Appends a range to a vector.
template <typename T, std::ranges::range Range> constexpr std::vector<T> operator+(const std::vector<T> &a, Range &&b)
{
    auto result = a;
    return result += std::forward<Range>(b);
}

/// Appends a range to a set.
template <typename T, std::ranges::range Range> constexpr std::set<T> operator+(const std::set<T> &a, Range &&b)
{
    auto result = a;
    return result += std::forward<Range>(b);
}
} // namespace euler
