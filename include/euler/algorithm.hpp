#pragma once

#include <map>
#include <random>
#include <set>

#include <boost/unordered/unordered_flat_map.hpp>

#include "counting_iterator.hpp"
#include "decls.hpp"
#include "it/base.hpp"

namespace euler
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
        *it++ = std::invoke(f, i);
    return result;
}

/// Returns the range [begin, end], possibly transformed by a function if given.
template <integral2 T, integral2 U, integral2 V, std::invocable<std::common_type_t<T, U, V>> Fun = std::identity>
constexpr auto range(T begin, U end, V step, Fun f = {})
{
    using Tp = std::common_type_t<T, U, V>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, Tp>>;

    Tp const b = begin;
    Tp const e = end;

    if (step == 0 || (step > 0 && e < b) || (step < 0 && e > b))
        return std::vector<FT>{};
    std::vector<FT> result((size_t)(step > 0 ? (e - b + step) / step : (b - e - step) / -step));
    auto it = result.begin();
    if (step > 0)
        for (Tp i = b; i <= e; i += step)
            *it++ = std::invoke(f, i);
    else
        for (Tp i = b; i >= e; i += step)
            *it++ = std::invoke(f, i);
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
        *it++ = std::invoke(f, i);
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
            *it++ = std::invoke(f, i);
    else
        for (int64_t i = Begin; i >= End; i += Step)
            *it++ = std::invoke(f, i);
    return result;
}

/// Unfolds a recursive sequence into a vector. The parameter `f` can take either the previous
/// element of the sequence or the pair (previous element, index) as parameter.
template <typename T, typename Fun> constexpr std::vector<T> unfold(size_t size, T seed, Fun f)
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
template <execution_policy Exec, std::ranges::range Range, typename Fun> auto mapv(Exec &&exec, Range &&r, Fun f)
{
    using T = std::ranges::range_value_t<Range>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
    std::vector<FT> result(r.size());
    std::transform(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), std::ranges::begin(result),
                   std::move(f));
    return result;
}

/// Transforms a range.
template <std::ranges::range Range, typename Fun> constexpr auto mapv(Range &&r, Fun f)
{
    using T = std::ranges::range_value_t<Range>;
    using FT = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
    std::vector<FT> result(r.size());
    std::transform(std::ranges::begin(r), std::ranges::end(r), std::ranges::begin(result), std::move(f));
    return result;
}

/// Reduces a range.
template <execution_policy Exec, std::ranges::range Range, typename T, typename BinaryOp, typename Fun = std::identity>
T reducev(Exec &&exec, Range &&r, T init, BinaryOp op, Fun f = {})
{
    return std::transform_reduce(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), std::move(init),
                                 std::move(op), std::move(f));
}

/// Reduces a range.
template <std::ranges::range Range, typename T, typename BinaryOp, typename Fun = std::identity>
constexpr T reducev(Range &&r, T init, BinaryOp op, Fun f = {})
{
    return std::transform_reduce(std::ranges::begin(r), std::ranges::end(r), std::move(init), std::move(op),
                                 std::move(f));
}

/// Reduces a range.
template <execution_policy Exec, integral2 T, integral2 U, typename Tp, typename BinaryOp, typename Fun = std::identity>
Tp reduceRange(Exec &&exec, T begin, U end, Tp init, BinaryOp op, Fun f = {})
{
    using V = std::common_type_t<T, U>;
    if (V(end) < V(begin))
        return init;
    return std::transform_reduce(std::forward<Exec>(exec), counting_iterator(V(begin)), counting_iterator(V(end + 1)),
                                 std::move(init), std::move(op), std::move(f));
}

/// Reduces a range.
template <integral2 T, integral2 U, typename Tp, typename BinaryOp, typename Fun = std::identity>
constexpr Tp reduceRange(T begin, U end, Tp init, BinaryOp op, Fun f = {})
{
    using V = std::common_type_t<T, U>;
    if (V(end) < V(begin))
        return init;
    return std::transform_reduce(counting_iterator(V(begin)), counting_iterator(V(end + 1)), std::move(init),
                                 std::move(op), std::move(f));
}

/// Sums a function over a range.
template <execution_policy Exec, std::ranges::range Range, typename Fun = std::identity>
auto sum(Exec &&exec, Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Exec>(exec), std::forward<Range>(r), T{}, std::plus{}, std::move(f));
}

/// Sums a function over a range.
template <std::ranges::range Range, typename Fun = std::identity> constexpr auto sum(Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Range>(r), T{}, std::plus{}, std::move(f));
}

/// Sums a function over a range of numbers.
template <execution_policy Exec, integral2 T, integral2 U, typename Fun = std::identity>
auto sum(Exec &&exec, T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end), Tp{}, std::plus{}, std::move(f));
}

/// Sums a function over a range of numbers.
template <integral2 T, integral2 U, typename Fun = std::identity> constexpr auto sum(T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::move(begin), std::move(end), Tp{}, std::plus{}, std::move(f));
}

/// Sums a function over a range of numbers using TBB.
template <integral2 T, integral2 U, typename Fun = std::identity> auto psum(T begin, U end, Fun f = {})
{
    using V = std::common_type_t<T, U>;
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return tbb::parallel_reduce(
        tbb::blocked_range<V>(begin, end + 1), Tp{},
        [&](tbb::blocked_range<V> r, Tp running_total) {
            for (auto i = r.begin(); i < r.end(); ++i)
                running_total += f(i);
            return running_total;
        },
        std::plus{});
}

template <size_t Threshold = 8192, integral2 T, integral2 U, typename Fun = std::identity>
auto sumMaybeParallel(T begin, U end, Fun f = {})
{
    if constexpr (Threshold == 0)
    {
        return sum(begin, end, f);
    }
    else
    {
        if (end - begin + 1 >= Threshold)
            return psum(begin, end, f);
        return sum(begin, end, f);
    }
}

/// Multiplies a function over a range.
template <execution_policy Exec, std::ranges::range Range, typename Fun = std::identity>
auto product(Exec &&exec, Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Exec>(exec), std::forward<Range>(r), T{1}, std::multiplies{}, std::move(f));
}

/// Multiplies a function over a range.
template <std::ranges::range Range, typename Fun = std::identity> constexpr auto product(Range &&r, Fun f = {})
{
    using T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Range>(r), T{1}, std::multiplies{}, std::move(f));
}

/// Multiplies a function over a range of numbers.
template <execution_policy Exec, integral2 T, integral2 U, typename Fun = std::identity>
auto product(Exec &&exec, T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end), Tp{1}, std::multiplies{},
                       std::move(f));
}

/// Multiplies a function over a range of numbers.
template <integral2 T, integral2 U, typename Fun = std::identity> constexpr auto product(T begin, U end, Fun f = {})
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, std::common_type_t<T, U>>>;
    return reduceRange(std::move(begin), std::move(end), Tp{1}, std::multiplies{}, std::move(f));
}

/// Returns the maximum value in a sequence according to a specified key selector function.
template <std::ranges::range Range, typename Key = std::identity>
constexpr auto maxBy(Range &&range, std::ranges::range_value_t<Range> init, Key key = {})
{
    return reducev(range, std::move(init), [&](auto &&a, auto &&b) {
        return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
    });
}

/// Returns the maximum value in a sequence according to a specified key selector function.
template <execution_policy Exec, std::ranges::range Range, typename Key = std::identity>
auto maxBy(Exec &&exec, Range &&range, std::ranges::range_value_t<Range> init, Key key = {})
{
    return reducev(std::forward<Exec>(exec), std::forward<Range>(range), std::move(init), [&](auto &&a, auto &&b) {
        return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
    });
}

/// Returns the maximum value in a numeric range according to a specified key selector function.
template <integral2 T, integral2 U, typename Tp, typename Key = std::identity>
constexpr auto maxRangeBy(T begin, U end, Tp init, Key key = {})
{
    return reduceRange(std::move(begin), std::move(end), std::move(init), [&](auto &&a, auto &&b) {
        return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
    });
}

/// Returns the maximum value in a numeric range according to a specified key selector function.
template <execution_policy Exec, integral2 T, integral2 U, typename Tp, typename Key = std::identity>
auto maxRangeBy(Exec &&exec, T begin, U end, Tp init, Key key = {})
{
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end), std::move(init),
                       [&](auto &&a, auto &&b) {
                           return key(a) < key(b) ? std::forward<decltype(b)>(b) : std::forward<decltype(a)>(a);
                       });
}

/// Returns a vector of partial sums of a function over a range.
template <execution_policy Exec, std::ranges::range Range, typename BinaryOp = std::plus<>,
          typename Fun = std::identity,
          typename T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>>
auto partialSum(Exec &&exec, Range &&r, T init = {}, BinaryOp op = {}, Fun f = {})
{
    std::vector<T> res(std::ranges::size(r));
    std::transform_inclusive_scan(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), std::begin(res),
                                  std::move(op), std::move(f), std::move(init));
    return res;
}

/// Returns a vector of partial sums of a function over a range.
template <std::ranges::range Range, typename BinaryOp = std::plus<>, typename Fun = std::identity,
          typename T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>>
constexpr auto partialSum(Range &&r, T init = {}, BinaryOp op = {}, Fun f = {})
{
    std::vector<T> res(std::ranges::size(r));
    std::transform_inclusive_scan(std::ranges::begin(r), std::ranges::end(r), std::begin(res), std::move(op),
                                  std::move(f), std::move(init));
    return res;
}

/// Takes partial sums of a function over a range in place.
template <execution_policy Exec, std::ranges::range Range, typename BinaryOp = std::plus<>,
          typename Fun = std::identity,
          typename T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>>
void partialSumInPlace(Exec &&exec, Range &&r, T init = {}, BinaryOp op = {}, Fun f = {})
{
    std::transform_inclusive_scan(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), std::begin(r),
                                  std::move(op), std::move(f), std::move(init));
}

/// Takes partial sums of a function over a range in place.
template <std::ranges::range Range, typename BinaryOp = std::plus<>, typename Fun = std::identity,
          typename T = std::remove_cvref_t<std::invoke_result_t<Fun, std::ranges::range_value_t<Range>>>>
constexpr void partialSumInPlace(Range &&r, T init = {}, BinaryOp op = {}, Fun f = {})
{
    std::transform_inclusive_scan(std::ranges::begin(r), std::ranges::end(r), std::begin(r), std::move(op),
                                  std::move(f), std::move(init));
}

/// Returns a vector of adjacent differences of a function over a range.
template <execution_policy Exec, std::ranges::range Range, typename BinaryOp = std::minus<>>
auto adjacentDifference(Exec &&exec, Range &&r, BinaryOp op = {})
{
    std::vector<std::ranges::range_value_t<Range>> res(std::ranges::size(r));
    std::adjacent_difference(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r), res.begin(),
                             std::move(op));
    return res;
}

/// Returns a vector of adjacent differences of a function over a range.
template <std::ranges::range Range, typename BinaryOp = std::minus<>>
constexpr auto adjacentDifference(Range &&r, BinaryOp op = {})
{
    std::vector<std::ranges::range_value_t<Range>> res(std::ranges::size(r));
    std::adjacent_difference(std::ranges::begin(r), std::ranges::end(r), res.begin(), std::move(op));
    return res;
}

/// Takes adjacent difference in place.
template <execution_policy Exec, std::ranges::range Range, typename BinaryOp = std::minus<>>
void adjacentDifferenceInPlace(Exec &&exec, Range &&r, BinaryOp op = {})
{
    std::adjacent_difference(std::forward<Exec>(exec), std::ranges::begin(r), std::ranges::end(r),
                             std::ranges::begin(r), std::move(op));
}

/// Takes adjacent difference in place.
template <std::ranges::range Range, typename BinaryOp = std::minus<>>
constexpr void adjacentDifferenceInPlace(Range &&r, BinaryOp op = {})
{
    std::adjacent_difference(std::ranges::begin(r), std::ranges::end(r), std::ranges::begin(r), std::move(op));
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

/// Returns the first integer `n` in `[low, high)` where `pred(n)` is true. Assumes that all false values precede all
/// true values. Returns `high` if no such integer exists.
template <typename U, std::predicate<U> Pred> constexpr U find_first(U low, U high, Pred pred)
{
    // Handle empty range or invalid input
    if (high <= low)
        return high;
    while (low < high)
    {
        U const mid = low + (high - low) / 2;
        if (pred(mid))
            high = mid;
        else
            low = mid + 1;
    }
    return low;
}

/// Returns the least integer `n` in `[low, high)` such that `f(n) ≥ target`, or `high` if no such integer exists.
/// Requires `f` to be nondecreasing.
template <typename T, typename U, std::invocable<U> Fun> constexpr U first_ge(Fun f, const T &target, U low, U high)
{
    return find_first(std::move(low), std::move(high), [&](const U &n) { return f(n) >= target; });
}

/// Returns the least integer `n` in `[low, high)` such that `f(n) > target`, or `high` if no such integer exists.
/// Requires `f` to be nondecreasing.
template <typename T, typename U, std::invocable<U> Fun> constexpr U first_gt(Fun f, const T &target, U low, U high)
{
    return find_first(std::move(low), std::move(high), [&](const U &n) { return f(n) > target; });
}

/// Returns the greatest integer `n` in `[low, high)` such that `f(n) ≤ target`, or `low - 1` if no such integer exists.
/// Requires `f` to be nondecreasing.
template <typename T, typename U, std::invocable<U> Fun> constexpr U last_le(Fun f, const T &target, U low, U high)
{
    return first_gt(std::move(f), std::move(target), std::move(low), std::move(high)) - 1;
}

/// Returns the greatest integer `n` in `[low, high)` such that `f(n) < target`, or `low - 1` if no such integer exists.
/// Requires `f` to be nondecreasing.
template <typename T, typename U, std::invocable<U> Fun> constexpr U last_lt(Fun f, const T &target, U low, U high)
{
    return first_ge(std::move(f), std::move(target), std::move(low), std::move(high)) - 1;
}

/// Returns the index of the first element in the sorted range `r` that is `≥ value`, or `size(r)` if no such element
/// exists.
template <std::ranges::random_access_range Range, typename T, typename Comp = std::ranges::less,
          typename Proj = std::identity>
[[nodiscard]] std::ranges::range_difference_t<Range> first_ge(Range &&r, const T &value, Comp comp = {}, Proj proj = {})
{
    return std::ranges::lower_bound(r, value, std::move(comp), std::move(proj)) - std::ranges::begin(r);
}

/// Returns the index of the first element in the sorted range `r` that is `> value`, or `size(r)` if no such element
/// exists.
template <std::ranges::random_access_range Range, typename T, typename Comp = std::ranges::less,
          typename Proj = std::identity>
[[nodiscard]] std::ranges::range_difference_t<Range> first_gt(Range &&r, const T &value, Comp comp = {}, Proj proj = {})
{
    return std::ranges::upper_bound(r, value, std::move(comp), std::move(proj)) - std::ranges::begin(r);
}

/// Returns the index of the last element in the sorted range `r` that is `≤ value`, or `-1` if no such element exists.
template <std::ranges::random_access_range Range, typename T, typename Comp = std::ranges::less,
          typename Proj = std::identity>
[[nodiscard]] std::ranges::range_difference_t<Range> last_le(Range &&r, const T &value, Comp comp = {}, Proj proj = {})
{
    return std::ranges::upper_bound(r, value, std::move(comp), std::move(proj)) - std::ranges::begin(r) - 1;
}

/// Returns the index of the last element in the sorted range `r` that is `< value`, or `-1` if no such element exists.
template <std::ranges::random_access_range Range, typename T, typename Comp = std::ranges::less,
          typename Proj = std::identity>
[[nodiscard]] std::ranges::range_difference_t<Range> last_lt(Range &&r, const T &value, Comp comp = {}, Proj proj = {})
{
    return std::ranges::lower_bound(r, value, std::move(comp), std::move(proj)) - std::ranges::begin(r) - 1;
}

/// Finds a root of the function `f`.
template <typename T, std::invocable<T> Fun> constexpr T findRoot(Fun f, T low, T high, size_t max_steps = 1000)
{
    auto const fl = f(low);
    auto const fh = f(high);
    if (fl == 0)
        return low;
    if (fh == 0)
        return high;
    if ((fl > 0) == (fh > 0))
        throw std::domain_error("Function has same sign at endpoints");
    if (fl > 0)
        std::swap(low, high);
    while (max_steps-- > 0)
    {
        T const mid = low + (high - low) / 2;
        auto const y = f(mid);
        if (low == mid || high == mid || y == 0)
            return mid;
        if (y < 0)
            low = mid;
        else
            high = mid;
    }
    return high;
}

/// Minimizes a unimodal function `f` on a real interval `[low, high]`.
/// @param max_steps Controls precision: error ~ (high - low) * (2/3)^iterations.
template <typename T, std::invocable<T> Fun> T ternarySearch(Fun f, T low, T high, size_t max_steps = 200)
{
    while (max_steps-- > 0)
    {
        T const m1 = low + (high - low) / 3;
        T const m2 = high - (high - low) / 3;
        if (m1 == m2)
            break;
        if (f(m2) < f(m1))
            low = m1; // minimum is in [m1, high]
        else
            high = m2; // minimum is in [low, m2]
    }
    return (low + high) / 2;
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
    std::vector<T> res;
    res.reserve(totalSize(v));
    for (const auto &e : std::forward<Range>(v))
        flatten(e, res);
    return res;
}

/// Same as `collections.Counter` in Python.
template <template <typename...> typename Map = std::map, std::ranges::range Range> auto counter(Range &&r)
{
    using T = std::ranges::range_value_t<Range>;
    Map<T, size_t> res;
    for (auto &&i : std::forward<Range>(r))
        ++res[i];
    return res;
}

/// Creates a histogram (frequency count) of elements in the range `v`.
template <std::ranges::range Range, std::integral T = std::ranges::range_value_t<Range>>
std::vector<T> histogram(const Range &r, T maxItem)
{
    std::vector<T> res(maxItem + 1);
    for (auto &&x : r)
        if (x >= 0 && x <= maxItem)
            ++res[x];
    return res;
}

/// Creates a histogram (frequency count) of elements in the range `v`.
template <execution_policy Exec, std::ranges::range Range, std::integral T = std::ranges::range_value_t<Range>>
std::vector<T> histogram(Exec && /*exec*/, Range &&r, T maxItem)
{
    if constexpr (std::same_as<std::decay_t<Exec>, std::execution::sequenced_policy> ||
                  std::same_as<std::decay_t<Exec>, std::execution::unsequenced_policy>)
        return histogram(std::forward<Range>(r), maxItem);
    else
        return tbb::parallel_reduce(
            tbb::blocked_range(std::ranges::begin(r), std::ranges::end(r)), std::vector<T>(maxItem + 1),
            [&](auto &&range, std::vector<T> &&h) {
                for (auto i = range.begin(); i != range.end(); ++i)
                    if (*i >= 0 && *i <= maxItem)
                        ++h[*i];
                return std::move(h);
            },
            [&](std::vector<T> a, const std::vector<T> &b) {
                for (T i = 0; i <= maxItem; ++i)
                    a[i] += b[i];
                return a;
            });
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
    period_result res{.period = 1, .preperiod = 0};
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
    using T = decltype(auto(fn(Z(0))));
    T total = sum(start, std::min(stop, preperiod - 1), [&](auto &&i) -> decltype(auto) { return fn(i); });
    T mid = 0;
    start = std::max(preperiod, start);
    Z const l = (stop - start + 1) % period + start - 1;
    for (Z i = start; i <= std::min(stop, start + period - 1); ++i)
    {
        auto y = fn(i);
        mid += y;
        if (i <= l)
            total += y;
    }
    return total + (stop - start + 1) / period * mid;
}

/// Convenience function that returns v, sorted and with duplicates removed.
template <std::ranges::range Range, typename Comp = std::ranges::less> auto sortedSet(Range &&r, Comp comp = {})
{
    std::vector v(std::ranges::begin(r), std::ranges::end(r));
    std::ranges::sort(v, comp);
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

/// Convenience function that returns v, sorted and with duplicates removed.
template <execution_policy Exec, std::ranges::range Range, typename Comp = std::ranges::less>
auto sortedSet(Exec &&exec, Range &&r, Comp comp = {})
{
    std::vector v(std::ranges::begin(r), std::ranges::end(r));
    std::sort(exec, v.begin(), v.end(), comp);
    v.erase(std::unique(std::forward<decltype(exec)>(exec), v.begin(), v.end()), v.end());
    return v;
}

/// Counting sort.
template <std::ranges::range Range, std::integral T = std::ranges::range_value_t<Range>>
void countSort(Range &&v, T maxItem)
{
    auto const h = histogram(v, maxItem);
    auto it = v.begin();
    for (size_t i = 0; i < h.size(); ++i)
        for (T k = 0; k < h[i]; ++k)
            *it++ = i;
}

/// Counting sort.
template <execution_policy Exec, std::ranges::range Range, std::integral T = std::ranges::range_value_t<Range>>
void countSort(Exec &&exec, Range &&v, T maxItem)
{
    auto const h = partialSum(histogram(std::forward<Exec>(exec), v, maxItem));
    it::range(0, h.size() - 1)(std::forward<Exec>(exec), [&](size_t i) {
        auto it = i == 0 ? v.begin() : v.begin() + h[i - 1];
        auto end = v.begin() + h[i];
        for (; it != end; ++it)
            *it = i;
    });
}

/// Convenience method to add a range to a vector.
template <typename T, std::ranges::range Range> std::vector<T> operator+(std::vector<T> v, Range &&r)
{
    v.append_range(r);
    return v;
}

/// Convenience method to add a range to a vector.
template <typename T, std::ranges::range Range> std::vector<T> &operator+=(std::vector<T> &v, Range &&r)
{
    v.append_range(r);
    return v;
}

/// Convenience method to repeat a vector n times.
template <typename T> std::vector<T> operator*(const std::vector<T> &v, size_t n)
{
    std::vector<T> result(n * v.size());
    for (size_t i = 0; i < n * v.size(); ++i)
        result[i] = v[i % v.size()];
    return result;
}

/// Convenience method to repeat a vector n times.
template <typename T> std::vector<T> operator*(const T &scalar, const std::vector<T> &a) { return a * scalar; }

// ==== Atomic helpers ====

/// Atomically replaces the current value with the result of std::min of the value and arg. That is, it performs atomic
/// minimum operation.
template <typename T> T fetch_min(std::atomic<T> &a, T b, std::memory_order m = std::memory_order_seq_cst) noexcept
{
    return __atomic_fetch_min((T *)&a, b, (int)m); // NOLINT
}

/// Atomically replaces the current value with the result of std::max of the value and arg. That is, it performs atomic
/// maximum operation.
template <typename T> T fetch_max(std::atomic<T> &a, T b, std::memory_order m = std::memory_order_seq_cst) noexcept
{
    return __atomic_fetch_max((T *)&a, b, (int)m); // NOLINT
}

/// Sums a function over the nodes of a binary tree. The `Node` type must have functions `left`, `right`, and `height`.
/// The `height` function is compared with `limit` to determine when to stop.
template <typename Node, typename T, typename Fun> auto sumBinaryTree(Node root, T limit, Fun f)
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, Node>>;
    if (root.height() > limit)
        return Tp{};
    std::vector<Node> st{root};
    Tp res{};
    while (!st.empty())
    {
        auto node = st.back();
        st.pop_back();
        for (; node.height() <= limit; node = *node.left())
        {
            res += f(node);
            auto const right = node.right();
            if (right && right->height() <= limit)
                st.push_back(*right);
            if (!node.left())
                break;
        }
    }
    return res;
}

/// Sums a function over the nodes of a binary tree in parallel. The `Node` type must have functions `left`, `right`,
/// and `height`. The `height` function is compared with `limit` to determine when to stop.
template <typename Node, typename T, typename U, typename Fun>
auto sumBinaryTree(Node root, T limit, U parallelThreshold, Fun f)
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, Node>>;
    if (root.height() > parallelThreshold)
        return sumBinaryTree(root, limit, f);

    Tp left_res{}, right_res{};
    auto const left = root.left(), right = root.right();
    bool const do_left = left && left->height() <= limit, do_right = right && right->height() <= limit;
    if (do_left && do_right)
        tbb::parallel_invoke([&] { right_res = sumBinaryTree(*right, limit, parallelThreshold, f); },
                             [&] { left_res = sumBinaryTree(*left, limit, parallelThreshold, f); });
    else if (do_left)
        left_res = sumBinaryTree(*left, limit, parallelThreshold, f);
    else if (do_right)
        right_res = sumBinaryTree(*right, limit, parallelThreshold, f);
    return f(root) + left_res + right_res;
}

/// Constructs a map from each element of `r` to its index in `r`.
template <std::ranges::random_access_range Range> [[nodiscard]] constexpr auto indexMap(Range &&r)
{
    using Key = std::ranges::range_value_t<Range>;
    using Index = std::ranges::range_size_t<Range>;
    auto const n = std::ranges::size(r);
    boost::unordered_flat_map<Key, Index> res;
    res.reserve(n);
    for (Index i = 0; i < n; ++i)
        res[r[i]] = i;
    return res;
}
} // namespace euler
