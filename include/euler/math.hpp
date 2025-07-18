#pragma once

#include "PF.hpp"
#include "modular_arithmetic.hpp"
#include <algorithm>
#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <numeric>

inline namespace euler
{
namespace detail
{
template <std::ranges::range Range> constexpr size_t countSquarefreeHelper(Range &&primes, size_t limit, size_t start)
{
    size_t result = 0;
    for (size_t i = start; i < primes.size(); ++i)
    {
        size_t const p = primes[i] * primes[i];
        size_t const l = limit / p;
        if (l == 0)
            break;
        result += l;
        if (l > p)
            result -= countSquarefreeHelper(primes, l, i + 1);
    }
    return result;
}
} // namespace detail

/// The factorials from 0 to 20, precisely those which fit in a 64-bit integer.
inline constexpr std::array<int64_t, 21> factorialsTo20{1,
                                                        1,
                                                        2,
                                                        6,
                                                        24,
                                                        120,
                                                        720,
                                                        5040,
                                                        40320,
                                                        362880,
                                                        3628800,
                                                        39916800,
                                                        479001600,
                                                        6227020800,
                                                        87178291200,
                                                        1307674368000,
                                                        20922789888000,
                                                        355687428096000,
                                                        6402373705728000,
                                                        121645100408832000,
                                                        2432902008176640000};

/// Calculates ⌊a / b⌋, and works for negative values too.
template <integral2 Ta, integral2 Tb> constexpr auto floorDiv(const Ta &a, const Tb &b)
{
    auto const d = boost::multiprecision::detail::evaluate_if_expression(a / b);
    return d * b == a ? d : d - ((a < 0) ^ (b < 0));
}

/// Calculates ⌈a / b⌉, and works for negative values too.
template <integral2 Ta, integral2 Tb> constexpr auto ceilDiv(const Ta &a, const Tb &b)
{
    return floorDiv(a + b - (b > 0 ? 1 : -1), b);
}

template <integral2 T> constexpr bool isSquare(const T &n)
{
    T a = isqrt(n);
    return a * a == n;
}

template <integral2 T> constexpr bool isSquare(const boost::rational<T> &r)
{
    return isSquare(r.numerator()) && isSquare(r.denominator());
}

/// Greatest common divisor of a range.
template <std::ranges::range Range> constexpr std::ranges::range_value_t<Range> gcd(Range &&r)
{
    using T = std::ranges::range_value_t<Range>;
    using std::gcd;
    T g = 0;
    for (auto &&x : r)
    {
        g = gcd(g, x);
        if (g == 1)
            break;
    }
    return g;
}

/// Greatest common divisor, specialized to an initializer list.
template <integral2 T> constexpr T gcd(std::initializer_list<T> l)
{
    return gcd<std::initializer_list<T>>(std::move(l));
}

/// Least common multiple of a range.
template <std::ranges::range Range> constexpr std::ranges::range_value_t<Range> lcm(Range &&r)
{
    using T = std::ranges::range_value_t<Range>;
    using std::lcm;
    T l = 1;
    for (auto &&x : r)
        l = lcm(l, x);
    return l;
}

/// Least common multiple, specialized to an initializer list.
template <integral2 T> constexpr T lcm(std::initializer_list<T> l)
{
    return lcm<std::initializer_list<T>>(std::move(l));
}

/// A replacement for integer floor division.
template <typename T, typename U> constexpr T fastDiv(T a, U b)
{
    T q = T((double)a / (double)b);
    T const prod = q * b;
    if (b > 0)
    {
        if (prod > a)
        {
            --q;
            while (q * b > a)
                --q;
        }
        else if (prod + b <= a)
        {
            ++q;
            while ((q + 1) * b <= a)
                ++q;
        }
    }
    else if (prod < a)
    {
        --q;
        while (q * b < a)
            --q;
    }
    else if (prod + b >= a)
    {
        ++q;
        while ((q + 1) * b >= a)
            ++q;
    }
    return q;
}

/// Wrapper for `mpz_divexact`.
inline mpz_int &divexact(mpz_int &dest, const mpz_int &a, const mpz_int &b)
{
    mpz_divexact(dest.backend().data(), a.backend().data(), b.backend().data());
    return dest;
}

/// Wrapper for `mpz_divexact_ui`.
inline mpz_int &divexact(mpz_int &dest, const mpz_int &a, uint32_t b)
{
    mpz_divexact_ui(dest.backend().data(), a.backend().data(), b);
    return dest;
}

/// Returns the square root of an integer, if it is a square. Otherwise, returns none.
template <integral2 T> constexpr std::optional<T> sqrtIfSquare(T n)
{
    auto s = isqrt(n);
    return s * s == n ? std::optional{s} : std::nullopt;
}

/// Returns the square root of a rational number, if it is a square. Otherwise, returns none.
template <integral2 T> constexpr std::optional<boost::rational<T>> sqrt(const boost::rational<T> &r)
{
    T num = isqrt(r.numerator());
    T denom = isqrt(r.denominator());
    if (num * num != r.numerator() || denom * denom != r.denominator())
        return std::nullopt;
    return {{num, denom}};
}

/// Returns a table of binomial coefficients of size `size`.
template <typename T = int64_t> constexpr std::vector<std::vector<T>> binomTable(size_t size)
{
    std::vector<std::vector<T>> table(size + 1);
    table[0] = {1};
    for (size_t n = 1; n <= size; ++n)
    {
        table[n] = std::vector<T>(n + 1);
        table[n][0] = T(1);
        table[n][n] = T(1);
        for (size_t r = 1; r < n; ++r)
            table[n][r] = table[n - 1][r - 1] + table[n - 1][r];
    }
    return table;
}

/// Returns a vector of binomial coefficients (n choose k) for k from 0 to `limit`.
template <typename T> std::vector<T> binomialVec(size_t n, size_t limit)
{
    std::vector<T> res(limit + 1);
    res[0] = 1;
    for (size_t i = 1; i <= limit; ++i)
        res[i] = res[i - 1] * (n - i + 1) / i;
    return res;
}

/// Returns a solution to two simultaneous congruences along with the lcm of the moduli. Requirements:
/// * `m, n > 0`.
/// * `lcm(m, n)` must be less than the maximum integer size.
/// There is a solution iff `gcd(m, n) | b - a`. If an invalid input is provided or there is no solution, returns
/// (`-1`, `-1`).
template <integral2 Ta, integral2 Tb, integral2 Tm, integral2 Tn> constexpr auto crtlcm(Ta a, Tb b, Tm m, Tn n)
{
    using T = decltype(auto(boost::multiprecision::detail::evaluate_if_expression(a * b * m * n)));
    if (m < 0 || n < 0)
        return std::pair<T, T>{-1, -1};
    auto const [g, s, _] = xgcd(m, n);
    T const diff = b - a;
    if (diff % g != 0)
        return std::pair<T, T>{-1, -1};
    T const l = m / g * n;
    return std::pair<T, T>{mod(a + modmul(modmul(diff / g, s, l), m, l), l), l};
}

/// Returns a solution to two or more simultaneous congruences along with the lcm of the moduli. Requirements:
/// * `m > 0` for all `m` ∈ `moduli`.
/// * `lcm(moduli)` must be less than the maximum integer size.
/// If an invalid input is provided or there is no solution, returns (`-1`, `-1`).
template <std::ranges::random_access_range Range1, std::ranges::random_access_range Range2>
constexpr auto crtlcm(Range1 &&remainders, Range2 &&moduli)
{
    using T = std::ranges::range_value_t<Range2>;

    return it::range(0UZ, std::min(remainders.size(), moduli.size()) - 1)
        .map([&](size_t i) { return std::pair<T, T>{remainders[i], moduli[i]}; })
        .reduce(std::pair<T, T>{0, 1},
                [&](auto &&a, auto &&b) { return crtlcm(a.first, b.first, a.second, b.second); });
}

/// Returns a solution to two simultaneous congruences. Requirements:
/// * `m, n > 0`.
/// * `lcm(m, n)` must be less than the maximum integer size.
/// There is a solution iff `gcd(m, n) | b - a`. If an invalid input is provided or there is no solution, returns `-1`.
template <integral2 T1, integral2 T2, integral2 T3, integral2 T4> constexpr auto crt(T1 a, T2 b, T3 m, T4 n)
{
    return crtlcm(a, b, m, n).first;
}

/// Returns a solution to two or more simultaneous congruences. Requirements:
/// * `m > 0` for all `m` ∈ `moduli`.
/// * `lcm(moduli)` must be less than the maximum integer size.
/// If an invalid input is provided or there is no solution, returns `-1`.
template <std::ranges::random_access_range Range1, std::ranges::random_access_range Range2>
constexpr auto crt(Range1 &&remainders, Range2 &&moduli)
{
    return crtlcm(std::forward<Range1>(remainders), std::forward<Range2>(moduli)).first;
}

/// Convenience method to return a list of powers of `a` from 0 to `n`.
template <typename T> constexpr std::vector<T> powers(T a, size_t n)
{
    std::vector<T> result(n + 1);
    T p = 1;
    for (auto &&x : result)
    {
        x = p;
        p *= a;
    }
    return result;
}

/// @brief Counts the number of integers in a certain residue class.
/// @param a The residue class.
/// @param modulus The modulus.
/// @param low The lower bound (inclusive).
/// @param high The upper bound (inclusive).
/// @return The number of elements in the residue class given by a.
template <integral2 Ta, integral2 Tm, integral2 Tl, integral2 Th>
constexpr auto countResidueClass(Ta a, Tm modulus, Tl low, Th high)
{
    // Alternative formula: floorDiv(high - a, modulus) + floorDiv(a - low, modulus) + 1;
    return floorDiv(high - a, modulus) - floorDiv(low - a - 1, modulus);
}

/// Returns `n!`.
template <typename T = int64_t> constexpr T factorial(int n)
{
    assert(n >= 0);
    // If T is std::integral, do assertions to fail if the result overflows.
    if constexpr (std::integral<T> && sizeof(T) <= 4)
        assert(n <= 12);
    else if constexpr (std::integral<T> && sizeof(T) <= 8)
        assert(n <= 20);
    else if constexpr (std::same_as<T, int128_t>)
        assert(n <= 33);
    else if constexpr (std::same_as<T, uint128_t>)
        assert(n <= 34);
    return product(1, n, [](int i) { return T(i); });
}

/// Returns `n!`.
template <> inline mpz_int factorial(int n)
{
    assert(n >= 0);
    mpz_int result;
    mpz_fac_ui((mpz_ptr)result.backend().data(), n);
    return result;
}

/// Returns `n!!`.
template <typename T = int64_t> constexpr T doubleFactorial(int n)
{
    assert(n >= 0);
    if constexpr (std::same_as<T, mpz_int>)
    {
        mpz_int result;
        mpz_2fac_ui((mpz_ptr)result.backend().data(), n);
        return result;
    }
    else
    {
        return product(0, floorDiv(n - 1, 2), [&](int i) { return T(n - 2 * i); });
    }
}

/// Calculates the binomial coefficient for given values of `n` and `k`.
template <typename T = int64_t> constexpr T binomial(int n, int r)
{
    if (r <= 0 || n == 0)
        return T(r == 0);
    if (n < 0)
        return pow(-1, r) * binomial<T>(r - n - 1, r);
    if (r > n)
        return 0;
    if constexpr (std::same_as<T, mpz_int>)
    {
        mpz_int result;
        mpz_bin_uiui((mpz_ptr)result.backend().data(), n, r);
        return result;
    }
    else
    {
        r = std::min(r, n - r); // Take advantage of symmetry
        T result{1};
        for (int i = 1; i <= r; ++i)
        {
            result *= n - r + i;
            result /= i;
        }
        return result;
    }
}

/// Calculates the binomial coefficient for given values of `n` and `k` using GMP.
inline mpz_int binomial(const mpz_int &n, int r)
{
    if (r <= 0 || n == 0)
        return {r == 0};
    if (n < 0)
        return pow(-1, r) * binomial(r - n - 1, r);
    if (r > n)
        return 0;
    mpz_int result;
    mpz_bin_ui((mpz_ptr)result.backend().data(), (mpz_srcptr)n.backend().data(), r);
    return result;
}

/// Calculates modular binomial coefficient using Lucas's theorem.
template <integral2 Tn, integral2 Tr, integral2 Tp> auto binomialMod(Tn n, Tr r, Tp p)
{
    using T = std::common_type_t<Tn, Tr, Tp>;
    T num = 1, denom = 1;
    while (n != 0 && r != 0)
    {
        auto const nmodp = n % p;
        auto rmodp = r % p;
        if (nmodp < rmodp)
            return T{};
        rmodp = std::min(rmodp, nmodp - rmodp);
        for (T i = 1; i <= rmodp; ++i)
        {
            num = num * (nmodp - rmodp + i) % p;
            denom = denom * i % p;
        }
        n /= p;
        r /= p;
    }
    return T(num * modInverse(denom, p) % p);
}

/// Calculates `Σ k ≥ 0, f(⌊n / (a * k + b)⌋)`.
template <execution_policy ExecutionPolicy, std::integral Tn, std::integral Ta = int, std::integral Tb = int,
          typename Fun = std::identity>
    requires(!std::integral<ExecutionPolicy>)
auto sumFloors(ExecutionPolicy &&exec, Tn n, Ta a = 1, Tb b = 1, Fun &&f = {})
{
    using T = std::common_type_t<Tn, Ta, Tb>;
    assert(a > 0 && b > 0);
    if constexpr (std::is_same_v<Fun, std::identity>)
        if (a == 1 && b == 1)
        {
            // Half as much computation in this branch.
            T sqrtn = isqrt(n);
            return 2 * sum(std::forward<ExecutionPolicy>(exec), T(1), sqrtn, [&](auto &&i) { return n / i; }) -
                   sqrtn * sqrtn;
        }
    auto K = std::max(T(0), std::min(n / b - 1, (T)isqrt(n)));
    auto res = sum(std::forward<ExecutionPolicy>(exec), T(1), K,
                   [&](auto &&k) { return f(k) * ((n - b * k) / (a * k) - (n - b * (k + 1)) / (a * k + a)); });
    if (n >= b * (K + 1))
        res += sum(std::forward<ExecutionPolicy>(exec), T(0), (n - b * (K + 1)) / (a * K + a),
                   [&](auto &&i) { return f(n / (a * i + b)); });
    return res;
}

/// Calculates `Σ k ≥ 0, f(⌊n / (a * k + b)⌋)`.
template <std::integral Tn, std::integral Ta = int, std::integral Tb = int, typename F = std::identity>
constexpr auto sumFloors(Tn n, Ta a = 1, Tb b = 1, F &&f = {})
{
    return sumFloors(std::execution::seq, n, a, b, std::forward<F>(f));
}

/// Calculates `∑ k ∈ [start, stop], f(⌊n / k⌋)`.
template <execution_policy Exec, std::integral T, typename Fun = std::identity>
T sumFloorsRange(Exec &&exec, T n, T start, T stop, Fun f = {})
{
    stop = std::min(n, stop);
    if (start > stop)
        return T(0);
    T c = isqrt(n);
    auto result = sum(std::forward<Exec>(exec), start, std::min(stop, c), [&](T k) { return f(n / k); });
    if (c < stop)
    {
        T a = n / stop;
        T b = n / std::max(c + 1, start);
        result += f(a) * (std::min(stop, n / a) - std::max(start - 1, n / (a + 1)));
        if (a != b)
            result += f(b) * (std::min(stop, n / b) - std::max(start - 1, n / (b + 1)));
        result += sum(std::forward<Exec>(exec), a + 1, b - 1, [&](T nk) { return f(nk) * (n / nk - n / (nk + 1)); });
    }
    return result;
}

/// Calculates `∑ k ∈ [start, stop], f(⌊n / k⌋)`.
template <std::integral T, typename Fun = std::identity> constexpr T sumFloorsRange(T n, T start, T stop, Fun f = {})
{
    stop = std::min(n, stop);
    if (start > stop)
        return T(0);
    T c = isqrt(n);
    auto result = sum(start, std::min(stop, c), [&](T k) { return f(n / k); });
    if (c < stop)
    {
        T a = n / stop;
        T b = n / std::max(c + 1, start);
        result += f(a) * (std::min(stop, n / a) - std::max(start - 1, n / (a + 1)));
        if (a != b)
            result += f(b) * (std::min(stop, n / b) - std::max(start - 1, n / (b + 1)));
        result += sum(a + 1, b - 1, [&](T nk) { return f(nk) * (n / nk - n / (nk + 1)); });
    }
    return result;
}

/// Number of integral `(x,y)` satisfying `x^2 + y^2 ≤ n`.
template <execution_policy ExecutionPolicy, std::integral T = int64_t>
auto latticePointsInDisk(ExecutionPolicy &&exec, const T &n)
{
    return 1 + 4 * (sumFloors(std::forward<ExecutionPolicy>(exec), n, T(4), T(1)) -
                    sumFloors(std::forward<ExecutionPolicy>(exec), n, T(4), T(3)));
}

/// Number of integral `(x,y)` satisfying `x^2 + y^2 ≤ n`.
auto latticePointsInDisk(const std::integral auto &n) { return latticePointsInDisk(std::execution::seq, n); }

/// Number of integral `(x,y)` satisfying `x^2 + y^2 ≤ n`. Somehow faster than hyperbola method.
auto latticePointsInDisk2(integral2 auto n)
{
    using T = decltype(n);
    T l = isqrt(n);
    T sum1;
    if (n < 600'000'000'000)
        sum1 = sum(T(1), l, [&](T x) -> T { return isqrt(n - x * x); });
    else
        sum1 = sum(std::execution::par, T(1), l, [&](T x) -> T { return isqrt(n - x * x); });
    return 1 + 4 * sum1 + 4 * l;
}

/// Counts squarefree numbers up to `n`. The mobius sieve needs to be filled up to `√n`.
template <execution_policy Exec, integral2 T> auto countSquarefree(Exec &&exec, T n, const std::vector<int8_t> &mobius)
{
    T sqrtn = isqrt(n);
    return sum(std::forward<Exec>(exec), 1, sqrtn,
               [&](integral2 auto i) { return mobius[i] == 0 ? 0 : mobius[i] * (n / (i * i)); });
}

/// Counts squarefree numbers up to `n`. The mobius sieve needs to be filled up to `√n`.
template <integral2 T> auto countSquarefree(T n, const std::vector<int8_t> &mobius)
{
    return n < 40'000'000 ? countSquarefree(std::execution::seq, n, mobius)
                          : countSquarefree(std::execution::par, n, mobius);
}

/// Counts squarefree numbers up to `n`. This is faster if you don't already have a mobius sieve.
template <integral2 T, std::ranges::range Range>
    requires(!std::same_as<std::ranges::range_value_t<Range>, int8_t>)
constexpr auto countSquarefree(T n, Range &&primes)
{
    return n - detail::countSquarefreeHelper(std::forward<Range>(primes), n, 0);
}

/// Calculates the nth polygonal number of a given type.
template <integral2 T> constexpr T polygonalNumber(int p, T n) { return (n * n * (p - 2) - n * (p - 4)) / 2; }

/// Sums `1 + 2 + ... + limit`, maximally avoiding overflow and also avoiding divisions in `T`.
template <typename T = int64_t> constexpr T sumId(size_t limit) { return T((limit + (limit & 1)) >> 1) * (limit | 1); }

/// Sums `1^2 + 2^2 + ... + limit^2`.
template <typename T = int64_t> constexpr T sumSquares(size_t limit)
{
    return sumId<T>(limit) * (2 * limit + 1) / 3;
    // Ugly code for maximally avoiding overflow and also avoiding divisions in `T`.
    // switch (limit % 6)
    // {
    // case 0:
    //     return T(limit / 6) * (limit + 1) * T(2 * limit + 1);
    // case 1:
    //     return T(limit) * T((limit + 1) / 2) * T((2 * limit + 1) / 3);
    // case 2:
    //     return T(limit / 2) * T((limit + 1) / 3) * T(2 * limit + 1);
    // case 3:
    //     return T(limit / 3) * T((limit + 1) / 2) * T(2 * limit + 1);
    // case 4:
    //     return T(limit / 2) * T(limit + 1) * T((2 * limit + 1) / 3);
    // case 5:
    //     return T(limit) * T((limit + 1) / 6) * T(2 * limit + 1);
    // default:
    //     std::unreachable();
    // }
}

/// Finds the largest `e` such that `b^e ≤ n`.
template <integral2 T, integral2 U> constexpr int floor_log(T n, U b)
{
    if (n < b)
        return 0;
    if (n < b * b)
        return 1;
    int e = 1;
    T x = b;
    for (; mulLeq(x, b, n); x *= b, ++e)
        ;
    return e;
}
} // namespace euler
