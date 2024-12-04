#pragma once

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
template <integral2 T> void doSpfSieveParallel(T start, T limit, T sqrtLimit, std::vector<T> &sieve)
{
    for (T i = start; i <= sqrtLimit; i += 2)
        if (sieve[i] == i)
            it::range(i * i, limit, i * 2)(std::execution::par, [&](T j) {
                if (sieve[j] == j)
                    sieve[j] = i;
            });
}

template <integral2 T> void doSpfSieveSequential(T start, T limit, T sqrtLimit, std::vector<T> &sieve)
{
    for (T i = start; i <= sqrtLimit; i += 2)
        if (sieve[i] == i)
            it::range(i * i, limit, i * 2)([&](T j) {
                if (sieve[j] == j)
                    sieve[j] = i;
            });
}

template <std::ranges::range Range> size_t countSquarefreeHelper(Range &&primes, size_t limit, size_t start)
{
    size_t result = 0;
    for (size_t i = start; i < primes.size(); ++i)
    {
        size_t p = primes[i] * primes[i];
        auto l = limit / p;
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
constexpr std::array<int64_t, 21> factorialsTo20{1LL,
                                                 1LL,
                                                 2LL,
                                                 6LL,
                                                 24LL,
                                                 120LL,
                                                 720LL,
                                                 5040LL,
                                                 40320LL,
                                                 362880LL,
                                                 3628800LL,
                                                 39916800LL,
                                                 479001600LL,
                                                 6227020800LL,
                                                 87178291200LL,
                                                 1307674368000LL,
                                                 20922789888000LL,
                                                 355687428096000LL,
                                                 6402373705728000LL,
                                                 121645100408832000LL,
                                                 2432902008176640000LL};

/// Factors `n!`.
template <integral2 T = int> constexpr std::vector<PrimePower<T>> factorFactorial(int n)
{
    return it::primes{2, n}.map([&](auto &&p) { return PrimePower<T>{p, factorialValuation(n, p)}; }).to();
}

/// Factors `binomial(n, r)`.
template <integral2 T = int> constexpr std::vector<PrimePower<T>> factorBinomial(int n, int r)
{
    return it::primes{2, n}
        .map([&](auto &&p) {
            auto e = factorialValuation(n, p) - factorialValuation(n - r, p) - factorialValuation(r, p);
            return PrimePower<T>{p, e};
        })
        .to();
}

/// Combines two prime factorizations using `op` on the exponents.
template <std::ranges::range Range, std::invocable<int, int> BinaryOp>
constexpr std::remove_cvref_t<Range> combinePF(const Range &r1, const std::ranges::range auto &r2, BinaryOp op)
{
    std::remove_cvref_t<Range> result;
    auto it1 = r1.begin();
    auto it2 = r2.begin();
    while (it1 != r1.end() && it2 != r2.end())
    {
        auto [p1, e1] = *it1;
        auto [p2, e2] = *it2;
        if (p1 == p2)
        {
            auto e = op(e1, e2);
            if (e != 0)
                result.emplace_back(p1, e);
            ++it1;
            ++it2;
        }
        else if (p1 < p2)
        {
            result.emplace_back(p1, op(e1, 0));
            ++it1;
        }
        else
        {
            result.emplace_back(p2, op(0, e2));
            ++it2;
        }
    }
    auto it3 = std::transform(it2, r2.end(), std::back_inserter(result),
                              [&](auto &&p) { return std::pair{p.first, op(0, p.second)}; });
    std::copy(it1, r1.end(), it3);
    return result;
}

template <std::ranges::range Range1, std::ranges::range Range2>
constexpr std::remove_cvref_t<Range1> mulPF(Range1 &&r1, Range2 &&r2)
{
    return combinePF(std::forward<Range1>(r1), std::forward<Range2>(r2), std::plus());
}

template <std::ranges::range Range1, std::ranges::range Range2>
constexpr std::remove_cvref_t<Range1> lcmPF(Range1 &&r1, Range2 &&r2)
{
    return combinePF(std::forward<Range1>(r1), std::forward<Range2>(r2), maximum());
}

template <std::ranges::range Range1, std::ranges::range Range2>
constexpr std::remove_cvref_t<Range1> gcdPF(Range1 &&r1, Range2 &&r2)
{
    return combinePF(std::forward<Range1>(r1), std::forward<Range2>(r2), minimum());
}

template <std::ranges::range Range1, std::ranges::range Range2>
constexpr std::remove_cvref_t<Range1> divPF(Range1 &&r1, Range2 &&r2)
{
    return combinePF(std::forward<Range1>(r1), std::forward<Range2>(r2), std::minus());
}

template <std::ranges::range Range> constexpr std::remove_cvref_t<Range> powPF(Range &&r, const int n)
{
    auto f2 = std::forward<Range>(r);
    for (auto &&x : f2)
        x.second *= n;
    return f2;
}

/// Returns the product of the prime factorization.
template <std::ranges::range Range> constexpr auto evalPF(Range &&r)
{
    return product(std::forward<Range>(r), [](auto &&t) { return pow(t.first, t.second); });
}

/// Returns a sieve containing smallest prime factors up to `limit`.
template <integral2 T> constexpr std::vector<T> spfSieve(T limit)
{
    std::vector<T> sieve(limit + 1);
    for (T i = 2; i <= limit; ++i)
        if (i % 2 == 0)
            sieve[i] = 2;
        else if (i % 3 == 0)
            sieve[i] = 3;
        else if (i % 5 == 0)
            sieve[i] = 5;
        else
            sieve[i] = i;

    T sqrtLimit = isqrt(limit);
    if (limit >= 2'000'000)
        detail::doSpfSieveParallel(T(7), T(limit), sqrtLimit, sieve);
    else
        detail::doSpfSieveSequential(T(7), T(limit), sqrtLimit, sieve);
    return sieve;
}

/// Sieve for divisor counts. This is faster than divisorCountSieve2 for limits over 2 million.
template <std::integral T> constexpr std::vector<T> divisorCountSieve(T limit)
{
    std::vector<T> sieve(limit + 1);
    for (T i = 1; i <= limit; ++i)
        sieve[i] = 1;
    for (T i = 2; i <= limit; ++i)
        if (sieve[i] == 1)
            for (T j = i; j <= limit; j += i)
                sieve[j] *= valuation(j, i) + 1;
    return sieve;
}

// Sieve for divisor counts. This is faster than divisorCountSieve for limits up to 2 million.
template <std::integral T> constexpr std::vector<T> divisorCountSieve2(T limit)
{
    std::vector<T> sieve(limit + 1);
    for (T i = 1; i <= limit; ++i)
        for (T j = i; j <= limit; j += i)
            ++sieve[j];
    return sieve;
}

/// Sieve for the σ₁ function, the divisor sum function.
template <std::integral T> constexpr std::vector<T> divisorSumSieve(T limit)
{
    std::vector<T> sieve(limit + 1);
    for (T i = 1; i <= limit; ++i)
        for (T j = i; j <= limit; j += i)
            sieve[j] += i;
    return sieve;
}

/// Sieve for the ω function, the number of distinct prime factors of a number.
constexpr std::vector<uint8_t> omegaSieve(size_t limit)
{
    std::vector<uint8_t> sieve(limit + 1);
    for (size_t i = 2; i <= limit; ++i)
        if (sieve[i] == 0)
            for (size_t j = i; j <= limit; j += i)
                ++sieve[j];
    return sieve;
}

/// Generates a sieve of the mobius function.
template <std::ranges::range Range>
constexpr std::vector<int8_t> mobiusSieve(const Range &spfs, size_t limit = std::numeric_limits<size_t>::max())
{
    limit = std::min(limit, spfs.size() - 1);
    std::vector<int8_t> μ(limit + 1);
    μ[1] = 1;
    for (size_t i = 2; i <= limit; ++i)
    {
        auto k = i / spfs[i];
        if (spfs[k] == spfs[i])
            μ[i] = 0;
        else
            μ[i] = -μ[k];
    }
    return μ;
}

/// Sieve for the Mobius function.
template <integral2 T> constexpr std::vector<int8_t> mobiusSieve(T limit)
{
    return mobiusSieve(spfSieve(limit), limit);
}

/// Generates a sieve of squarefree numbers up to a given limit.
constexpr std::vector<bool> squarefreeSieve(size_t limit)
{
    std::vector sieve(limit + 1, true);
    sieve[0] = false;
    size_t s = isqrt(limit);
    for (size_t i = 2; i <= s; ++i)
        if (sieve[i])
            for (size_t j = i * i; j <= limit; j += i * i)
                sieve[j] = false;
    return sieve;
}

/// Sieve for the totient function.
template <std::integral T> constexpr std::vector<T> totientSieve(T limit)
{
    std::vector<T> sieve(limit + 1);
    sieve[1] = 1;
    for (T i = 2; i <= limit; ++i)
        if (sieve[i] == 0)
        {
            sieve[i] = i - 1;
            for (T j = 2; i * j <= limit; ++j)
            {
                if (sieve[j] == 0)
                    continue;

                T q = j;
                T f = i - 1;
                while (q % i == 0)
                {
                    f *= i;
                    q /= i;
                }
                sieve[i * j] = f * sieve[q];
            }
        }
    return sieve;
}

/// Returns (s, t) where d = s²t and t is squarefree.
template <integral2 T, std::ranges::range Range> constexpr std::pair<T, T> sqfreeDecompose(T num, Range &&factorization)
{
    T s = 1;
    for (auto &&[p, e] : std::forward<Range>(factorization))
    {
        auto q = pow(p, e / 2);
        s *= q;
        num /= q * q;
    }
    return {s, num};
}

/// Returns (s, t) where d = s²t and t is squarefree.
template <integral2 T> constexpr std::pair<T, T> sqfreeDecompose(T num) { return sqfreeDecompose(num, factor(num)); }

template <integral2 T = int64_t> constexpr std::vector<std::vector<T>> binomTable(int size, const T &modulus = 0)
{
    std::vector<std::vector<T>> table(size + 1);
    table[0] = {1};
    for (int n = 1; n <= size; ++n)
    {
        table[n] = std::vector<T>(n + 1);
        table[n][0] = 1;
        table[n][n] = 1;
        for (int r = 1; r < n; ++r)
        {
            table[n][r] = table[n - 1][r - 1] + table[n - 1][r];
            if (modulus > 0)
                table[n][r] %= modulus;
        }
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
    auto &&[g, x, _] = extendedEuclidean(m, n);
    T diff = b - a;
    if (diff % g != 0)
        return std::pair<T, T>{-1, -1};
    T l = m / g * n;
    return std::pair{mod(T(a + modmul(modmul(diff / g, x, l), m, l)), l), l};
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
        .reduce(std::execution::unseq, std::pair<T, T>{0, 1},
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

// Function to find smallest primitive root of p
template <integral2 T, std::ranges::range Range = std::vector<int>>
constexpr std::common_type_t<int64_t, T> primitiveRoot(const T &p, Range &&spfs = {})
{
    using Tp = std::common_type_t<int64_t, T>;
    if (p == 2)
        return Tp(1);
    auto phi = p - 1;
    auto pf = factor(phi, std::forward<Range>(spfs));
    for (Tp r = 2; r <= phi; ++r)
    {
        bool flag = false;
        for (const auto &[q, e] : pf)
        {
            if (powm(r, phi / q, p) == 1)
            {
                flag = true;
                break;
            }
        }

        if (!flag)
            return r;
    }
    assert(false && "primitiveRoot: Should not reach here (maybe p wasn't prime).");
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

/// Computes `∑ (ab ≤ n), g(a) * h(b)`.
/// Also known as `∑ (k ≤ n), (g * h)(k)`.
/// Also known as `∑ (k ≤ n), g(k) * H(n/k) = ∑ (k ≤ n), h(k) * G(n/k)`.
/// Requirements:
/// * `G` and `H` are the summatory functions of `g` and `h`
/// * Need to be able to evaluate `g(k)`, `h(k)` for `k ≤ √n` and `G(m)` and `H(m)` for `m ≥ √n`.
template <execution_policy Exec, std::integral T, std::invocable<T> Fun1, std::invocable<T> SummatoryFun1,
          std::invocable<T> Fun2, std::invocable<T> SummatoryFun2>
auto sumConvolution(Exec &&exec, T n, Fun1 g, SummatoryFun1 G, Fun2 h, SummatoryFun2 H)
{
    T sqrtn = isqrt(n);
    return sum(std::forward<Exec>(exec), T(1), sqrtn, [&](auto i) { return g(i) * H(n / i); }) +
           sum(std::forward<Exec>(exec), T(1), sqrtn, [&](auto i) { return h(i) * G(n / i); }) - G(sqrtn) * H(sqrtn);
}

/// Computes `∑ (ab ≤ n), g(a) * h(b)`.
/// Also known as `∑ (k ≤ n), (g * h)(k)`.
/// Also known as `∑ (k ≤ n), g(k) * H(n/k) = ∑ (k ≤ n), h(k) * G(n/k)`.
/// Requirements:
/// * `G` and `H` are the summatory functions of `g` and `h`
/// * Need to be able to evaluate `g(k)`, `h(k)` for `k ≤ √n` and `G(m)` and `H(m)` for `m ≥ √n`.
template <std::integral T, std::invocable<T> Fun1, std::invocable<T> SummatoryFun1, std::invocable<T> Fun2,
          std::invocable<T> SummatoryFun2>
auto sumConvolution(T n, Fun1 g, SummatoryFun1 G, Fun2 h, SummatoryFun2 H)
{
    return n < 40'000'000 ? sumConvolution(std::execution::seq, n, g, G, h, H)
                          : sumConvolution(std::execution::par, n, g, G, h, H);
}

/// Calculates `Σ k ֫≥ 0, f(⌊n / (a * k + b)⌋)`.
template <execution_policy ExecutionPolicy, std::integral Tn, std::integral Ta = int, std::integral Tb = int,
          typename Fun = std::identity>
    requires(!std::integral<ExecutionPolicy>)
constexpr auto sumFloors(ExecutionPolicy &&exec, Tn n, Ta a = 1, Tb b = 1, Fun &&f = {})
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

/// Calculates `Σ k ֫≥ 0, f(⌊n / (a * k + b)⌋)`.
template <std::integral Tn, std::integral Ta = int, std::integral Tb = int, typename F = std::identity>
constexpr auto sumFloors(Tn n, Ta a = 1, Tb b = 1, F &&f = {})
{
    return sumFloors(std::execution::seq, n, a, b, std::forward<F>(f));
}

/// Calculates `∑ k ∈ [start, stop], f(⌊n / k⌋)`.
template <execution_policy Exec, std::integral T, typename Fun = std::identity>
constexpr T sumFloorsRange(Exec &&exec, T n, T start, T stop, Fun f = {})
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
constexpr auto latticePointsInDisk(ExecutionPolicy &&exec, const T &n)
{
    return 1 + 4 * (sumFloors(std::forward<ExecutionPolicy>(exec), n, T(4), T(1)) -
                    sumFloors(std::forward<ExecutionPolicy>(exec), n, T(4), T(3)));
}

/// Number of integral `(x,y)` satisfying `x^2 + y^2 ≤ n`.
constexpr auto latticePointsInDisk(const std::integral auto &n) { return latticePointsInDisk(std::execution::seq, n); }

/// Number of integral `(x,y)` satisfying `x^2 + y^2 ≤ n`. Somehow faster than hyperbola method.
constexpr auto latticePointsInDisk2(integral2 auto n)
{
    using T = decltype(n);
    T l = isqrt(n);
    T sum1;
    if (n < 600'000'000'000LL)
        sum1 = sum(T(1), l, [&](T x) -> T { return isqrt(n - x * x); });
    else
        sum1 = sum(std::execution::par, T(1), l, [&](T x) -> T { return isqrt(n - x * x); });
    return 1 + 4 * sum1 + 4 * l;
}

/// Counts squarefree numbers up to `n`. The mobius sieve needs to be filled up to `√n`.
template <execution_policy Exec, integral2 T> auto countSquarefree(Exec &&exec, T n, const std::vector<int8_t> &mobius)
{
    T sqrtn = isqrt(n);
    return sum(std::forward<Exec>(exec), 1LL, sqrtn,
               [&](integral2 auto i) { return mobius[i] == 0 ? 0 : mobius[i] * (n / (i * i)); });
}

/// Counts squarefree numbers up to `n`. The mobius sieve needs to be filled up to `√n`.
template <integral2 T> constexpr auto countSquarefree(T n, const std::vector<int8_t> &mobius)
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

/// An efficient sparse array whose keys are floors ⌊N/i⌋ for i >= 1.
/// Usage: `floors_array(n)` or `floors_array(n, s)`.
template <integral2 T = int64_t> class floors_array
{
  public:
    /// s should be 1 less than the size of the small values array.
    constexpr floors_array(size_t n, size_t s) : _n(n), _up(s + 1), _down(n / s) {}
    constexpr explicit floors_array(size_t n) : floors_array(n, isqrt(n)) {}

    /// Stores `∑ i ∈ [1..k], (μ * h)(i) = ∑ i ∈ [1..k], μ(i) * H(⌊k/i⌋)` for each `k` in the floors array. Here,
    /// `H` is the summatory function of `h`, and `hr` is the range of precomputed values of `h`.
    template <std::invocable<int64_t> SummatoryFun, std::ranges::random_access_range Range = std::array<T, 1>>
    [[nodiscard]] static floors_array mu(size_t n, SummatoryFun H, Range &&hr = {0})
    {
        size_t s = std::min(n / 2, std::max(isqrt(n), hr.size() - 1));
        s = n / ceilDiv(n, s); // Align on boundary
        size_t const ms = std::min(hr.size(), s + 1);
        floors_array fa(n, s);
        std::inclusive_scan(hr.begin() + 1, hr.begin() + ms, fa._up.begin() + 1, std::plus{}, T(0));
        for (size_t k = ms; k <= s; ++k)
            fa._up[k] = H(k) - sumFloors(k, 1, 2, [&](auto &&j) { return fa._up[j]; });
        for (size_t nk = fa._down.size() - 1; nk != 0; --nk)
        {
            auto k = (int64_t)(fa._n / nk);
            fa._down[nk] = H(k) - sumFloors(k, 1, 2, [&](auto &&j) { return fa[j]; });
        }
        return fa;
    }

    /// Stores `∑ i ∈ [1..k], φ(i)` for each `k` in the floors array. For best performance, precompute a totient sieve
    /// up to `0.4 * n^(2/3)` for this.
    template <std::ranges::random_access_range Range = std::array<T, 1>>
    [[nodiscard]] static floors_array sumTotient(size_t n, Range &&totients = {0})
    {
        return mu(n, [&](auto &&x) -> T { return T(x) * (x + 1) / 2; }, std::forward<Range>(totients));
    }

    constexpr T &operator[](size_t i) { return i < _up.size() ? _up[i] : _down[_n / i]; }
    constexpr const T &operator[](size_t i) const { return i < _up.size() ? _up[i] : _down[_n / i]; }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] constexpr size_t n() const { return _n; }
    /// The transition point between up and down. What was passed as the s parameter during construction.
    [[nodiscard]] constexpr size_t transitionPoint() const { return _up.size() - 1; }
    /// The up vector.
    [[nodiscard]] constexpr const std::vector<T> &up() const { return _up; }
    /// The down vector.
    [[nodiscard]] constexpr const std::vector<T> &down() const { return _down; }

    /// Enumerates keys of this floors array in ascending order. Breaks if `f` returns
    /// `it::result_break`.
    template <typename Fun> constexpr it::result_t ascending(Fun f) const
    {
        for (size_t i = 1; i < _up.size(); ++i)
            if (!it::callbackResult(f, i))
                return it::result_break;
        for (size_t i = _down.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, _n / i))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates keys of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <typename Fun> constexpr it::result_t descending(Fun f) const
    {
        for (size_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, _n / i))
                return it::result_break;
        for (size_t i = _up.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, i))
                return it::result_break;
        return it::result_continue;
    }

  private:
    size_t _n;
    std::vector<T> _up;
    std::vector<T> _down;
};

/// Calculates the sum of a completely multiplicative function applied to prime numbers in a given range, and returns
/// the computed map.
///
/// @return The map storing sum ∑ f(p) for p ∈ [2, n), p prime for all n of the form floor(limit/k).
template <std::invocable<int64_t> Fun, std::invocable<int64_t> SummatoryFun>
constexpr auto sumPrimeRangeLucyHedgehogExt(int64_t limit, Fun f, SummatoryFun F)
{
    using FT = std::remove_cvref_t<std::invoke_result_t<SummatoryFun, int64_t>>;
    int64_t r = isqrt(limit);
    auto fa = floors_array<FT>{static_cast<size_t>(limit)};
    fa.ascending([&](size_t i) { fa[i] = F(i) - F(1); });
    for (size_t p = 2; p <= (size_t)r; ++p)
    {
        if (fa[p] <= fa[p - 1])
            continue;
        // p is prime at this point.
        fa.descending([&](size_t i) {
            if (i < p * p)
                return it::result_break;
            fa[i] -= f(p) * (fa[i / p] - fa[p - 1]);
            return it::result_continue;
        });
    }
    return fa;
}

/// Calculates `∑ (p ∈ [2, limit] and p prime), f(p)`. Note that you must pass in the summatory function `F` of `f`, not
/// `f` itself.
template <std::invocable<int64_t> Fun, std::invocable<int64_t> SummatoryFun>
constexpr auto sumPrimeRangeLucyHedgehog(int64_t limit, Fun f, SummatoryFun F)
{
    return sumPrimeRangeLucyHedgehogExt(limit, std::move(f), std::move(F))[limit];
}

/// Calculates `∑ (p ∈ [2, limit] and p prime), p`.
template <typename T = int64_t> constexpr T sumPrimeRangeLucyHedgehog(int64_t limit)
{
    return sumPrimeRangeLucyHedgehog(limit, std::identity{}, [](auto &&n) { return T(n) * (n + 1) / 2; });
}
} // namespace euler
