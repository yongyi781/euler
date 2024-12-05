#pragma once

#include "algorithm.hpp"
#include "concepts.hpp"
#include "it/digits.hpp"
#include <boost/unordered/unordered_flat_map.hpp>
#include <ranges>

inline namespace euler
{
constexpr int64_t primePiUpperBound(int64_t n)
{
    if (n <= 20)
        return 8;
    int64_t result = (int64_t)((double)n / log(n) * (1.0 + 1.0 / log(n) + 2.51 / log(n) / log(n)));
    if (n <= 355990)
        return result + 5;
    return result;
}

/**
 * Enumerates the digits of a number.
 *
 * @tparam TBase the type of the base.
 * @param num The number.
 * @param base The base.
 * @param callback The callback function to invoke for each digit. Return false to break out of the loop.
 *
 * @return true if the loop terminates to completion, false otherwise
 */
template <integral2 TBase = int>
constexpr bool enumDigits(integral2 auto num, TBase base, std::invocable<TBase> auto callback)
{
    while (num > 0)
    {
        if (!invokeTrueIfVoid(callback, TBase(num % base)))
            return false;
        num /= base;
    }
    return true;
}

/// @brief Adds two numbers modulo another number in place.
/// @param a A number reference.
/// @param b A number.
constexpr void addm(integral2 auto &a, const integral2 auto &b, const integral2 auto &modulus)
{
    a = mod(a + b, modulus);
}

/// @brief Subtracts two numbers modulo another number in place.
/// @param a A number reference.
/// @param b A number.
constexpr void subm(integral2 auto &a, const integral2 auto &b, const integral2 auto &modulus)
{
    addm(a, modulus - b, modulus);
}

/// @brief Multiplies two numbers modulo another number in place.
/// @param a A number reference.
/// @param b A number.
/// @param modulus The modulus.
constexpr void mulm(integral2 auto &a, const integral2 auto &b, const integral2 auto &modulus)
{
    a = mod(a * b, modulus);
}

template <typename It, integral2 T>
constexpr bool _enumSmoothNumbersHelper(const T &limit, It it, It end, T current, std::invocable<T> auto callback)
{
    while (true)
    {
        auto limit2 = limit / current;
        if (!invokeTrueIfVoid(callback, current))
            return false;
        for (It it2 = std::next(it); it2 != end && *it2 <= limit2; ++it2)
            if (!_enumSmoothNumbersHelper(limit, it2, end, T(current * *it2), callback))
                return false;
        if (*it > limit2)
            break;
        current *= *it;
    }
    return true;
}

/**
 * Enumerates smooth numbers up to a given limit using a callback function.
 *
 * @tparam Range The type of the range storing prime numbers.
 * @tparam T The type of the limit.
 * @param primes The collection of prime numbers. It must be sorted.
 * @param limit The limit up to which to enumerate smooth numbers.
 * @param callback The callback function to be called for each smooth number.
 */
template <std::ranges::range Range, integral2 T>
bool enumSmoothNumbers(Range &&primes, T limit, std::invocable<T> auto callback)
{
    if (!invokeTrueIfVoid(callback, T(1)))
        return false;
    for (auto it = primes.begin(); it != primes.end() && *it <= limit; ++it)
        if (!_enumSmoothNumbersHelper(limit, it, primes.end(), T(*it), callback))
            return false;
    return true;
}

/**
 * Returns a vector of smooth numbers up to a given limit.
 *
 * @tparam Range The type of the range storing prime numbers.
 * @tparam T The type of the limit.
 * @param primes The collection of prime numbers. It must be sorted.
 * @param limit The limit up to which to enumerate smooth numbers.
 * @return A vector of smooth numbers up to the limit.
 */
template <std::ranges::range Range, integral2 T> constexpr std::vector<T> smoothNumbers(Range &&primes, T limit)
{
    std::vector<T> result;
    enumSmoothNumbers(std::forward<Range>(primes), limit, [&](auto &&n) { result.push_back(n); });
    return result;
}

template <typename It, integral2 T>
constexpr bool _enumSquarefreeSmoothNumbersHelper(const T &limit, It it, It end, T current,
                                                  std::invocable<T> auto callback)
{
    invokeTrueIfVoid(callback, current);
    auto limit2 = limit / current;
    for (It it2 = std::next(it); it2 != end && *it2 <= limit2; ++it2)
        if (!_enumSquarefreeSmoothNumbersHelper(limit, it2, end, T(current * *it2), callback))
            return false;
    return true;
}

/**
 * Enumerates squarefree smooth numbers up to a given limit using a callback function.
 *
 * @tparam Range The type of the range storing prime numbers.
 * @tparam T The type of the limit.
 * @param primes The collection of prime numbers. It must be sorted.
 * @param limit The limit up to which to enumerate squarefree smooth numbers.
 * @param callback The callback function to be called for each squarefree smooth number.
 */
template <std::ranges::range Range, integral2 T>
constexpr bool enumSquarefreeSmoothNumbers(Range &&primes, T limit, std::invocable<T> auto callback)
{
    callback(T(1));
    for (auto it = primes.begin(); it != primes.end() && *it <= limit; ++it)
        if (!_enumSquarefreeSmoothNumbersHelper(limit, it, primes.end(), T(*it), callback))
            return false;
    return true;
}

/**
 * Returns a vector of squarefree smooth numbers up to a given limit.
 *
 * @tparam Range The type of the range storing prime numbers.
 * @tparam T The type of the limit.
 * @param primes The collection of prime numbers. It must be sorted.
 * @param limit The limit up to which to enumerate squarefree smooth numbers.
 * @return A vector of squarefree smooth numbers up to the limit.
 */
template <std::ranges::range Range, integral2 T>
constexpr std::vector<T> squarefreeSmoothNumbers(Range &&primes, T limit)
{
    std::vector<T> result;
    enumSquarefreeSmoothNumbers(std::forward<Range>(primes), limit, [&](auto &&n) { result.push_back(n); });
    return result;
}

template <typename It, typename T>
constexpr bool _enumDivisorsHelper(It it, It end, const T &current, std::invocable<T> auto callback)
{
    if (!invokeTrueIfVoid(callback, current))
        return false;
    for (; it != end; ++it)
    {
        auto &&[p, e] = *it;
        T n = current;
        for (int i = 1; i <= e; ++i)
        {
            n *= p;
            if (!_enumDivisorsHelper(std::next(it), end, n, callback))
                return false;
        }
    }
    return true;
}

template <std::ranges::range Range>
constexpr bool enumDivisors(Range &&factorization,
                            std::invocable<typename std::ranges::range_value_t<Range>::first_type> auto callback)
{
    using T = std::ranges::range_value_t<Range>::first_type;
    return _enumDivisorsHelper(factorization.begin(), factorization.end(), T(1), callback);
}

template <std::ranges::range Range,
          std::invocable<typename std::ranges::range_value_t<Range>::first_type> UnaryOp = std::identity>
constexpr auto sumDivisors(Range &&factorization, UnaryOp f = std::identity())
{
    using T = std::ranges::range_value_t<Range>::first_type;
    if constexpr (std::is_same_v<UnaryOp, std::identity>)
    {
        return product(std::forward<Range>(factorization),
                       [](auto &&pe) { return (pow(pe.first, pe.second + 1) - 1) / (pe.first - 1); });
    }
    else
    {
        T result = 0;
        enumDivisors(std::forward<Range>(factorization), [&](auto &&d) { result += f(d); });
        return result;
    }
}

template <std::ranges::range Range = std::vector<int>>
constexpr auto sumDivisors(const integral2 auto &num, Range &&spfs = {})
{
    return reduceFactor(
        num, (decltype(num))1, std::multiplies(),
        [](auto &&pe) {
            if (pe.second <= 2)
                return sum(0, pe.second, [&](int i) { return pow(pe.first, i); });
            return (pow(pe.first, pe.second + 1) - 1) / (pe.first - 1);
        },
        std::forward<Range>(spfs));
}

template <std::ranges::range Range> constexpr auto countDivisors(Range &&factorization)
{
    return product(std::forward<Range>(factorization), [](auto &&pe) { return pe.second + 1; });
}

template <integral2 T, std::ranges::range Range = std::vector<int>>
constexpr auto countDivisors(const T &num, Range &&spfs = {})
{
    return reduceFactor(
        num, T(1), std::multiplies(), [](auto &&pe) { return T(pe.second + 1); }, std::forward<Range>(spfs));
}

/** Alternative (old) implementation of divisors. Faster in some cases. */
template <std::ranges::range Range> constexpr auto divisors2(Range &&factorization)
{
    using T = std::ranges::range_value_t<Range>::first_type;
    if (factorization.empty())
        return std::vector<T>{1};

    return cartesianFold(std::forward<Range>(factorization) |
                             std::views::transform([](auto &&pe) { return powers(pe.first, pe.second); }),
                         std::multiplies());
}

template <typename T = int, integral2 TBase = int, std::invocable<TBase> U = std::identity,
          std::invocable<T, std::remove_reference_t<std::invoke_result_t<U, TBase>>> B = std::plus<>>
constexpr T reduceDigits(const integral2 auto &num, TBase base, T identity, B binaryOp, U f = std::identity())
{
    T total = identity;
    it::digits(num, base)([&](TBase d) { total = binaryOp(total, f(d)); });
    return total;
}

/** Folds each tuple in the Cartesian product of a pair of lists using a binary operation. */
template <std::ranges::range Range1, std::ranges::range Range2,
          std::invocable<std::ranges::range_value_t<Range1>, std::ranges::range_value_t<Range2>> BinaryOp>
constexpr auto cartesianFold(Range1 &&r1, Range2 &&r2, BinaryOp op)
{
    using T = std::remove_cvref_t<
        std::invoke_result_t<BinaryOp, std::ranges::range_value_t<Range1>, std::ranges::range_value_t<Range2>>>;
    size_t size = r1.size() * r2.size();
    std::vector<T> result(size);
    auto it = result.begin();
    for (auto &&x : r1)
        for (auto &&y : r2)
            *it++ = T(op(x, y));
    return result;
}

/** Folds each tuple in the Cartesian product of a list of lists using a binary operation. */
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<std::ranges::range_value_t<Range>>,
                                                   std::ranges::range_value_t<std::ranges::range_value_t<Range>>>
                                        BinaryOp>
constexpr auto cartesianFold(Range &&ranges, BinaryOp op)
{
    assert(!ranges.empty());
    auto it = ranges.begin();
    auto result = *it++;
    for (; it != ranges.end(); ++it)
        result = cartesianFold(result, *it, op);
    return result;
}

template <std::ranges::range Range, typename It>
constexpr bool _enumCombinationsHelper(Range &combination, It it, It end, decltype(combination.begin()) cit,
                                       std::invocable<Range> auto callback)
{
    if (cit == combination.end())
        return invokeTrueIfVoid(callback, combination);
    while (it != end)
    {
        *cit = *it;
        if (!_enumCombinationsHelper(combination, ++it, std::next(end, 1), std::next(cit), callback))
            return false;
    }
    return true;
}

/// Enumerates all possible combinations of k elements from a range.
template <std::ranges::range Range>
constexpr bool enumCombinations(Range &&v, size_t k,
                                std::invocable<std::vector<std::ranges::range_value_t<Range>>> auto callback)
{
    if (k > v.size())
        return true; // True because iteration completed successfully on zero elements.

    std::vector<std::ranges::range_value_t<Range>> combination(k, 0);
    return _enumCombinationsHelper(combination, v.begin(), std::prev(v.end(), k - 1), combination.begin(), callback);
}

template <std::ranges::range Range, typename It>
constexpr bool _enumCombinationsWithReplacementHelper(Range &combination, It it, It end,
                                                      decltype(combination.begin()) cit,
                                                      std::invocable<Range> auto callback)
{
    if (cit == combination.end())
        return invokeTrueIfVoid(callback, combination);
    while (it != end)
    {
        *cit = *it;
        if (!_enumCombinationsWithReplacementHelper(combination, it++, end, std::next(cit), callback))
            return false;
    }
    return true;
}

/** Enumerates all possible combinations of k elements from a range, with replacement. */
template <std::ranges::range Range>
constexpr bool
enumCombinationsWithReplacement(Range &&v, size_t k,
                                std::invocable<std::vector<std::ranges::range_value_t<Range>>> auto callback)
{
    std::vector<std::ranges::range_value_t<Range>> combination(k, 0);
    return _enumCombinationsWithReplacementHelper(combination, v.begin(), v.end(), combination.begin(), callback);
}

/** Enumerates all possible permutations of k elements from a range. The range must be sorted. */
template <std::ranges::range Range>
constexpr bool enumPermutations(Range &&v, size_t k,
                                std::invocable<std::vector<std::ranges::range_value_t<Range>>> auto callback)
{
    // Pass by value is intentional here.
    return enumCombinations(v, k, [&](auto comb) {
        do
        {
            if (!invokeTrueIfVoid(callback, comb))
                return false;
        } while (std::next_permutation(comb.begin(), comb.end()));
        return true;
    });
}

/** The `k` fold self-Cartesian product of a range. Also known as permutations with replacement. */
template <std::ranges::range Range>
constexpr bool enumCartesianPower(Range &&v, size_t k,
                                  std::invocable<std::vector<std::ranges::range_value_t<Range>>> auto callback)
{
    // Pass by value is intentional here.
    return enumCombinationsWithReplacement(v, k, [&](auto comb) {
        do
        {
            if (!invokeTrueIfVoid(callback, comb))
                return false;
        } while (std::next_permutation(comb.begin(), comb.end()));
        return true;
    });
}

template <std::ranges::range Range, typename It>
constexpr bool _enumCartesianPowerSortedHelper(Range &&v, size_t k,
                                               std::vector<std::ranges::range_value_t<Range>> &combination, It cit,
                                               std::invocable<Range> auto callback)
{
    if (cit == combination.end())
        return invokeTrueIfVoid(callback, combination);
    for (auto &&x : v)
    {
        *cit = x;
        if (!_enumCartesianPowerSortedHelper(v, k, combination, std::next(cit), callback))
            return false;
    }
    return true;
}

/** The `k` fold self-Cartesian product of a range. This one is faster than non-sorted for smaller values of `k`. */
template <std::ranges::range Range>
constexpr bool enumCartesianPowerSorted(Range &&v, size_t k,
                                        std::invocable<std::vector<std::ranges::range_value_t<Range>>> auto callback)
{
    std::vector<std::ranges::range_value_t<Range>> combination(k);
    return _enumCartesianPowerSortedHelper(v, k, combination, combination.begin(), callback);
}

template <std::ranges::range Range, typename It>
constexpr bool _enumCartesianProductHelper(Range &combination, It it, It end, decltype(combination.begin()) cit,
                                           std::invocable<Range> auto callback)
{
    if (cit == combination.end())
        return invokeTrueIfVoid(callback, combination);
    for (auto &&x : *it)
    {
        *cit = x;
        if (!_enumCartesianProductHelper(combination, std::next(it), end, std::next(cit), callback))
            return false;
    }
    return true;
}

/** Cartesian product of a list of lists. */
template <std::ranges::range Range>
    requires std::ranges::range<std::ranges::range_value_t<Range>>
constexpr bool enumCartesianProduct(
    Range &&vs,
    std::invocable<std::vector<std::ranges::range_value_t<std::ranges::range_value_t<Range>>>> auto callback)
{
    std::vector<std::ranges::range_value_t<std::ranges::range_value_t<Range>>> combination(vs.size());
    return _enumCartesianProductHelper(combination, vs.begin(), vs.end(), combination.begin(), callback);
}

template <std::ranges::range Range>
    requires std::ranges::range<std::ranges::range_value_t<Range>>
constexpr auto cartesianProduct(Range &&vs)
{
    std::vector<std::vector<std::ranges::range_value_t<std::ranges::range_value_t<Range>>>> result;
    enumCartesianProduct(vs, [&](auto &&comb) { result.push_back(comb); });
    return result;
}

template <std::ranges::range Range, typename It>
constexpr bool _enumPowersetHelper(Range &combination, It it, It end, std::invocable<Range> auto callback)
{
    if (!invokeTrueIfVoid(callback, combination))
        return false;
    while (it != end)
    {
        combination.push_back(*it);
        if (!_enumPowersetHelper(combination, ++it, end, callback))
            return false;
        combination.pop_back();
    }
    return true;
}

/** Enumerates all subsets in a range. */
template <std::ranges::range Range>
constexpr bool enumPowerset(Range &&v, std::invocable<std::vector<std::ranges::range_value_t<Range>>> auto callback)
{
    std::vector<std::ranges::range_value_t<Range>> combination;
    combination.reserve(v.size());
    return _enumPowersetHelper(combination, v.begin(), v.end(), callback);
}

/**
 * Returns the powerset of a range.
 *
 * @tparam Range The type of the range.
 * @param v The input vector.
 * @return The powerset of the input vector.
 */
template <std::ranges::range Range> constexpr auto powerset(Range &&v)
{
    std::vector<std::vector<std::ranges::range_value_t<Range>>> result(pow(2uz, v.size()));
    auto it = result.begin();
    enumPowerset(v, [&](auto &&comb) { *it++ = comb; });
    return result;
}

constexpr bool _enumTuplesWithSumHelper(std::ranges::range auto &t, size_t index, int64_t total,
                                        std::invocable<std::vector<int64_t>> auto callback)
{
    for (int64_t i = 1; i <= total; ++i)
    {
        t[index] = i;
        t.back() = total - i;
        if (!invokeTrueIfVoid(callback, t))
            return false;
        if (total - i > 0)
            for (size_t index2 = t.size() - 2; index2 > index; index2--)
                if (!_enumTuplesWithSumHelper(t, index2, total - i, callback))
                    return false;
        t[index] = 0;
    }
    return true;
}

/** Enumerate nonnegative tuples with a given sum. Takes about 1-3.5 ns per tuple depending on length. */
constexpr bool enumTuplesWithSum(size_t length, int64_t total, std::invocable<std::vector<int64_t>> auto callback)
{
    std::vector t(length, (int64_t)0);
    t.back() = total;
    if (!invokeTrueIfVoid(callback, t))
        return false;
    for (size_t index = t.size() - 2; index != (size_t)(-1); index--)
        if (!_enumTuplesWithSumHelper(t, index, total, callback))
            return false;
    return true;
}

constexpr bool _enumPartitionsHelper(std::ranges::range auto &partition, int64_t n,
                                     std::invocable<std::vector<int64_t>> auto callback)
{
    partition.push_back(n);
    if (!invokeTrueIfVoid(callback, partition))
        return false;
    partition.pop_back();

    int64_t lb = partition.size() == 0 ? 1 : partition.back();
    for (int64_t i = n / 2; i >= lb; i--)
    {
        partition.push_back(i);
        if (!_enumPartitionsHelper(partition, n - i, callback))
            return false;
        partition.pop_back();
    }
    return true;
}

/** Enumerates all partitions of `n`. Takes about 4.8 ns per partition generated. */
constexpr bool enumPartitions(int64_t n, std::invocable<std::vector<int64_t>> auto callback)
{
    std::vector<int64_t> partition;
    partition.reserve(n);
    return _enumPartitionsHelper(partition, n, callback);
}

template <integral2 T, std::ranges::range Range = std::vector<int>>
constexpr bool enumFactor(const T &num, std::invocable<PrimePower<T>> auto callback, Range &&spfs = {})
{
    assert(num != 0 && "0 does not have a factorization");
    T n = num;
    if (n < 0)
    {
        if (!invokeTrueIfVoid(callback, PrimePower<T>{-1, 1}))
            return false;
        n = -n;
    }
    T start = 2;
    while (n > 1)
    {
        T p = !spfs.empty() && n < (int64_t)spfs.size() ? spfs[(size_t)n] : smallestPrimeFactor(n, start);
        if (!invokeTrueIfVoid(callback, PrimePower<T>{p, valuationDivide(n, p)}))
            return false;
        if (p == 2)
            start = 3;
        else
            start = p + 2;
    }
    return true;
}

template <integral2 T, std::ranges::range Range = std::vector<int>, std::invocable<PrimePower<T>> U,
          std::invocable<std::remove_reference_t<std::invoke_result_t<U, PrimePower<T>>>,
                         std::remove_reference_t<std::invoke_result_t<U, PrimePower<T>>>>
              B>
constexpr auto reduceFactor(const T &num, std::remove_cvref_t<std::invoke_result_t<U, PrimePower<T>>> identity,
                            B binaryOp, U f, Range &&spfs = {})
{
    auto total = identity;
    enumFactor(num, [&](auto &&pe) { total = binaryOp(total, f(pe)); }, std::forward<Range>(spfs));
    return total;
}

/// Returns the product of the elements in a range modulo a modulus.
template <execution_policy Exec, std::ranges::range Range,
          std::invocable<std::ranges::range_value_t<Range>> UnaryOp = std::identity>
constexpr auto modproduct(Exec &&exec, Range &&v, integral2 auto modulus, UnaryOp f = std::identity())
{
    using T = std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<decltype(exec)>(exec), std::forward<Range>(v), T{1}, mod_multiplies(modulus), f);
}

/// Returns the product of the elements in a range modulo a modulus.
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> UnaryOp = std::identity>
constexpr auto modproduct(Range &&v, integral2 auto modulus, UnaryOp f = std::identity())
{
    using T = std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Range>(v), T{1}, mod_multiplies(modulus), f);
}

/// Returns the product of the elements in a range modulo a modulus.
template <execution_policy Exec, integral2 TBegin, integral2 TEnd,
          std::invocable<std::common_type_t<TBegin, TEnd>> UnaryOp = std::identity>
constexpr auto modproduct(Exec &&exec, TBegin begin, TEnd end, integral2 auto modulus, UnaryOp f = std::identity())
{
    return reduceRange(std::forward<decltype(exec)>(exec), begin, end,
                       std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::common_type_t<TBegin, TEnd>>>{1},
                       mod_multiplies(modulus), f);
}

/// Returns the product of the elements in a range modulo a modulus.
template <integral2 TBegin, integral2 TEnd, std::invocable<std::common_type_t<TBegin, TEnd>> UnaryOp = std::identity>
constexpr auto modproduct(TBegin begin, TEnd end, integral2 auto modulus, UnaryOp f = std::identity())
{
    return reduceRange(begin, end,
                       std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::common_type_t<TBegin, TEnd>>>{1},
                       mod_multiplies(modulus), f);
}

/// Sums a function over a range, modulo a modulus.
template <execution_policy Exec, std::ranges::range Range,
          std::invocable<std::ranges::range_value_t<Range>> UnaryOp = std::identity>
constexpr auto modsum(Exec &&exec, Range &&v, integral2 auto modulus, UnaryOp f = std::identity())
{
    using T = std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Exec>(exec), std::forward<Range>(v), T{}, mod_plus(modulus), f);
}

/// Sums a function over a range, modulo a modulus.
template <std::ranges::range Range, std::invocable<std::ranges::range_value_t<Range>> UnaryOp = std::identity>
constexpr auto modsum(Range &&v, integral2 auto modulus, UnaryOp f = std::identity())
{
    using T = std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::ranges::range_value_t<Range>>>;
    return reducev(std::forward<Range>(v), T{}, mod_plus(modulus), f);
}

/// Sums a function over a range of values modulo a modulus.
template <execution_policy Exec, integral2 TBegin, integral2 TEnd,
          std::invocable<std::common_type_t<TBegin, TEnd>> UnaryOp = std::identity>
constexpr auto modsum(Exec &&exec, TBegin begin, TEnd end, integral2 auto modulus, UnaryOp f = std::identity())
{
    return reduceRange(std::forward<Exec>(exec), std::move(begin), std::move(end),
                       std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::common_type_t<TBegin, TEnd>>>{},
                       mod_plus(modulus), f);
}

/// Sums a function over a range of values modulo a modulus.
template <integral2 TBegin, integral2 TEnd, std::invocable<std::common_type_t<TBegin, TEnd>> UnaryOp = std::identity>
constexpr auto modsum(TBegin begin, TEnd end, integral2 auto modulus, UnaryOp f = std::identity())
{
    return reduceRange(begin, end,
                       std::remove_cvref_t<std::invoke_result_t<UnaryOp, std::common_type_t<TBegin, TEnd>>>{},
                       mod_plus(modulus), f);
}

auto _inverseConvolveSumHelper(std::ranges::range auto &&m, int64_t n, std::invocable<int64_t> auto H,
                               std::invocable<int64_t> auto g, std::invocable<int64_t> auto G)
{
    int64_t sqrtn = isqrt(n);
    if (sqrtn > 200)
        return H(n) - sum(std::execution::par, 2LL, sqrtn, [&](int64_t i) { return g(i) * m[n / i]; }) -
               sum(std::execution::par, 1, n / sqrtn - 1,
                   [&](int64_t d) -> int64_t { return (G(n / d) - G(n / (d + 1))) * m[d]; });
    return H(n) - sum(2LL, sqrtn, [&](int64_t i) { return g(i) * m[n / i]; }) -
           sum(1, n / sqrtn - 1, [&](int64_t d) -> int64_t { return (G(n / d) - G(n / (d + 1))) * m[d]; });
}

// floors_array is a lot faster here.

/// Returns the map containing values of sum((h * g^-1)(k) for k <= n), for n of the form limit/i.
/// H is the summatory of h.
auto inverseConvolveSumTable(int64_t limit, std::invocable<int64_t> auto H, std::invocable<int64_t> auto g,
                             std::invocable<int64_t> auto G)
{
    using TResult = std::remove_reference_t<std::invoke_result_t<decltype(H), int64_t>>;
    int64_t sqrtn = isqrt(limit);
    boost::unordered_flat_map<int64_t, TResult> result;
    for (int64_t i = 1; i <= sqrtn; ++i)
        result[i] = _inverseConvolveSumHelper(result, i, H, g, G);
    for (int64_t i = limit / sqrtn - 1; i > 0; i--)
        result[limit / i] = _inverseConvolveSumHelper(result, limit / i, H, g, G);
    return result;
}

/// Returns the map containing values of sum(mu(k)H(n/k) for k <= n), for n of the form limit/i.
template <integral2 T>
auto mobiusSumTable(int64_t limit, const std::vector<T> &smallValues, std::invocable<int64_t> auto H)
{
    using TResult = std::remove_reference_t<std::invoke_result_t<decltype(H), int64_t>>;
    int64_t cutoff = std::min(std::max((int64_t)5e6, isqrt(limit)), (int64_t)smallValues.size() - 1);
    boost::unordered_flat_map<int64_t, TResult> result;
    for (int64_t i = 1; i <= cutoff; ++i)
        result[i] = result[i - 1] + smallValues[i];
    for (int64_t i = limit / cutoff - 1; i > 0; i--)
        result[limit / i] = _inverseConvolveSumHelper(result, limit / i, H, [](auto) { return T(1); }, std::identity());
    return result;
}

/// Returns the map containing values of mertens(n/k) for all k <= n.
inline auto mertensTable(int64_t limit, const std::vector<int8_t> &mobius)
{
    return mobiusSumTable(limit, mobius, [](auto) { return (int64_t)1; });
}

auto _mobiusModsumHelper(std::ranges::range auto &m, int64_t n, integral2 auto modulus, std::invocable<int64_t> auto f)
{
    int64_t sqrtn = isqrt(n);
    if (sqrtn > 200)
        return mod(f(n) - modsum(std::execution::par, 2, sqrtn, modulus, [&](int64_t i) { return m[n / i]; }) -
                       modsum(std::execution::par, 1, n / sqrtn - 1, modulus,
                              [&](int64_t d) -> int64_t { return (n / d - n / (d + 1)) * m[d]; }),
                   modulus);
    return mod(
        f(n) - modsum(2, sqrtn, modulus, [&](int64_t i) { return m[n / i]; }) -
            modsum(1, n / sqrtn - 1, modulus, [&](int64_t d) -> int64_t { return (n / d - n / (d + 1)) * m[d]; }),
        modulus);
}

/// Returns the map containing values of sum(mu(k)f(n/k) for k <= n), for n of the form limit/i.
template <integral2 T>
auto mobiusModsumTable(int64_t limit, const std::vector<T> &table, integral2 auto modulus,
                       std::invocable<int64_t> auto f)
{
    using TResult = std::remove_reference_t<std::invoke_result_t<decltype(f), int64_t>>;
    int64_t cutoff = std::min({(int64_t)5e6, limit, (int64_t)table.size() - 1});
    boost::unordered_flat_map<int64_t, TResult> result;
    for (int64_t i = 1; i <= cutoff; ++i)
        result[i] = (result[i - 1] + table[i]) % modulus;
    for (int64_t i = limit / cutoff - 1; i > 0; i--)
        result[limit / i] = _mobiusModsumHelper(result, limit / i, modulus, f) % modulus;
    return result;
}

template <integral2 T>
auto reducePythagoreanTriples(auto &&exec, int64_t limit, const T &identity, auto binaryOp, auto f)
{
    int64_t hu = isqrt(limit);
    return transform_reduce(std::forward<decltype(exec)>(exec), counting_iterator((int64_t)1),
                            counting_iterator(hu + 1), identity, binaryOp, [&](auto &&u) {
                                T result = identity;
                                for (T v = (u % 2) + 1; v <= u - 1; v += 2)
                                {
                                    if (std::gcd(u, v) != 1)
                                        continue;
                                    auto a = u * u - v * v;
                                    auto b = 2 * u * v;
                                    auto c = u * u + v * v;
                                    if (c > limit)
                                        break;
                                    auto maxK = limit / c;
                                    for (auto k = 1; k <= maxK; ++k)
                                        result = binaryOp(result, f(k * a, k * b, k * c));
                                }
                                return result;
                            });
}
} // namespace euler
