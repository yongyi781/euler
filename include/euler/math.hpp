#pragma once

#include "PF.hpp"
#include "decls.hpp"
#include "it/primes.hpp"
#include "modular_arithmetic.hpp"
#include "prime.hpp"
#include <boost/unordered/unordered_flat_map.hpp>
#include <numeric>
#include <primesieve.hpp>

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
constexpr std::array<int64_t, 21> factorialsTo20{1,
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

/// Factors `n!`.
template <integral2 T = int> constexpr PF<T> factorFactorial(int n)
{
    return it::primes{2, n}.map([&](auto &&p) { return PrimePower<T>{p, factorialValuation(n, p)}; }).to();
}

/// Factors `binomial(n, r)`.
template <integral2 T = int> constexpr PF<T> factorBinomial(int n, int r)
{
    return it::primes{2, n}
        .map([&](auto &&p) {
            auto e = factorialValuation(n, p) - factorialValuation(n - r, p) - factorialValuation(r, p);
            return PrimePower<T>{p, e};
        })
        .to();
}

/// Sieve for the totient function.
template <std::integral T> constexpr std::vector<T> totientSieve(T limit)
{
    std::vector<T> phi(limit + 1);
    std::vector<T> primes;
    phi[1] = 1;
    primes.reserve(limit / std::max(1.0, log(limit)));

    for (T i = 2; i <= limit; i++)
    {
        if (phi[i] == 0)
        {
            phi[i] = i - 1;
            primes.push_back(i);
        }
        for (T p : primes)
        {
            if (!mulLeq(i, p, limit))
                break;
            if (i % p == 0)
            {
                phi[i * p] = phi[i] * p;
                break;
            }
            phi[i * p] = phi[i] * (p - 1);
        }
    }
    return phi;
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

/// Returns a solution to two simultaneous congruences along with the lcm of the moduli. Requirements:
/// * `m, n > 0`.
/// * `lcm(m, n)` must be less than the maximum integer size.
/// There is a solution iff `gcd(m, n) | b - a`. If an invalid input is provided or there is no solution, returns
/// (`-1`, `-1`).
template <integral2 Ta, integral2 Tb, integral2 Tm, integral2 Tn> constexpr auto crtlcm(Ta a, Tb b, Tm m, Tn n)
{
    using T = std::common_type_t<Ta, Tb, Tm, Tn>;

    if (m < 0 || n < 0)
        return std::pair<T, T>{-1, -1};
    auto &&[g, x, _] = xgcd(m, n);
    T diff = b - a;
    if (diff % g != 0)
        return std::pair<T, T>{-1, -1};
    T l = m / g * n;
    return std::pair{mod(T(a + modmul(modmul(T(diff / g), x, l), m, l)), l), l};
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

template <typename T> constexpr std::vector<T> powers(T a, int n)
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
constexpr auto countResidueClass(integral2 auto a, integral2 auto modulus, integral2 auto low, integral2 auto high)
{
    // Alternative formula: floorDiv(high - a, modulus) + floorDiv(a - low, modulus) + 1;
    return floorDiv(high - a, modulus) - floorDiv(low - a - 1, modulus);
}

inline mpz_int factorial(int n)
{
    assert(n >= 0);
    mpz_int result;
    mpz_fac_ui((mpz_ptr)result.backend().data(), n);
    return result;
}

/// Calculates the binomial coefficient for given values of `n` and `k` using GMP.
inline mpz_int binomial(int n, int r)
{
    if (r == 0)
        return 1;
    if (n < 0 || r < 0 || r > n)
        return 0;
    mpz_int result;
    mpz_bin_uiui((mpz_ptr)result.backend().data(), n, r);
    return result;
}

/// Calculates the binomial coefficient for given values of `n` and `k` using GMP.
inline mpz_int binomial(const mpz_int &n, int r)
{
    if (r == 0)
        return 1;
    if (n < 0 || r < 0 || r > n)
        return 0;
    mpz_int result;
    mpz_bin_ui((mpz_ptr)result.backend().data(), (mpz_srcptr)n.backend().data(), r);
    return result;
}

/// Calculates the binomial coefficient for given values of n and k using a naive iterative approach.
template <integral2 T, integral2 U> constexpr auto binomial2(T n, U r)
{
    using Tp = std::common_type_t<T, U>;
    if (r > n)
        return Tp(0);
    if (r == 0 || r == n)
        return Tp(1);
    if (r > n - r)
        r = U(n - r); // Take advantage of symmetry
    Tp result{1};
    for (U i = 1; i <= r; ++i)
    {
        result *= n - r + i;
        result /= i;
    }
    return result;
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

/// Generate the B + 1 numbers (B(1) = 1/2) up to `k`.
template <typename T> std::vector<T> bernoulliPlus(size_t k)
{
    std::vector<T> A(k + 1), B(k + 1);
    for (size_t m = 0; m <= k; m++)
    {
        A[m] = T(1) / T(m + 1);
        for (int j = m; j >= 1; j--)
            A[j - 1] = j * (A[j - 1] - A[j]);
        if (m == 1 || m % 2 == 0)
            B[m] = A[0];
    }
    return B;
}

/// Sums `1 + 2 + ... + limit`, maximally avoiding overflow and also avoiding divisions in `T`.
template <typename T = int64_t> T sumId(size_t limit)
{
    return limit % 2 == 0 ? T(limit / 2) * (limit + 1) : T(limit) * ((limit + 1) / 2);
}

/// Sums `1^2 + 2^2 + ... + limit^2`, maximally avoiding overflow and also avoiding divisions in `T`.
template <typename T = int64_t> T sumSquares(size_t limit)
{
    assert(limit >= 0);
    switch (limit % 6)
    {
    case 0:
        return T(limit / 6) * (limit + 1) * T(2 * limit + 1);
    case 1:
        return T(limit) * T((limit + 1) / 2) * T((2 * limit + 1) / 3);
    case 2:
        return T(limit / 2) * T((limit + 1) / 3) * T(2 * limit + 1);
    case 3:
        return T(limit / 3) * T((limit + 1) / 2) * T(2 * limit + 1);
    case 4:
        return T(limit / 2) * T(limit + 1) * T((2 * limit + 1) / 3);
    case 5:
        return T(limit) * T((limit + 1) / 6) * T(2 * limit + 1);
    default:
        std::unreachable();
    }
}
} // namespace euler
