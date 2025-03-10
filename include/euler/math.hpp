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
template <std::integral T> constexpr std::vector<T> spfSieve(T limit)
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

// Space-optimized structure for smallest prime factors (SPF) up to n.
// It stores SPF only for odd numbers. For any even number (>2), the SPF is 2.
template <std::integral T> struct SPF
{
    // Only need to store integer half the size of the input!
    using H = std::make_unsigned_t<half_integer_t<T>>;
    // spfOdd[i] holds the smallest prime factor for number (2*i + 1).
    // Index 0 corresponds to 1 (unused), index 1 to 3, index 2 to 5, etc.
    std::vector<H> spfOdd;
    std::vector<H> smallPrimes; // Odd primes up to sqrt(n).

    SPF() = default;
    SPF(T n) : spfOdd((n + 1) / 2, 0), smallPrimes(primeRange(H(3), H(isqrt(n))))
    {
        T const m = (n + 1) / 2; // covers numbers 1,3,5,... up to n
        // We ignore index 0 (number 1) and process indices [1, m).
        // Partition these indices into segments.
        int const segSize = 32768; // you can adjust this block size as needed
        std::vector<std::pair<T, T>> segments;
        for (T start = 1; start < m; start += segSize)
            segments.emplace_back(start, std::min(m, start + segSize));

        // Process each segment in parallel.
        std::for_each(std::execution::par, segments.begin(), segments.end(), [&](const std::pair<T, T> &seg) {
            T const L = seg.first;
            T const R = seg.second; // R is not inclusive
            // The segment covers odd numbers from:
            //   segStart = 2*L + 1  up to  segEnd = 2*R + 1 (exclusive)
            T const segStart = 2 * L + 1;
            T const segEnd = 2 * R + 1;
            // For each small prime p, mark its odd multiples in this segment.
            for (H const p : smallPrimes)
            {
                // Compute the first number to mark in this segment.
                T startVal = p * p;
                // If p*p is not less than the upper bound of the segment, nothing to mark.
                if (startVal >= segEnd)
                    break;
                if (startVal < segStart)
                {
                    // Increase startVal in steps of 2*p until it's >= segStart.
                    T const diff = segStart - startVal;
                    T const k = (diff + 2 * p - 1) / (2 * p); // ceil(diff / (2*p))
                    startVal += k * 2 * p;
                }
                // Mark every odd multiple of p in [startVal, segEnd)
                for (T x = startVal; x < segEnd; x += 2 * p)
                {
                    T const idx = (x - 1) / 2;
                    // Only update if not already marked (ensuring the smallest prime factor remains)
                    if (spfOdd[idx] == 0)
                        spfOdd[idx] = p;
                }
            }
        });
    }

    /// Accessor: returns the smallest prime factor for any x (1 <= x <= n).
    [[nodiscard]] T get(T n) const
    {
        if (n < 2)
            return 0;
        if (n % 2 == 0)
            return 2;
        return spfOdd[n / 2] == 0 ? n : spfOdd[n / 2];
    }

    /// Returns whether the SPF sieve is empty.
    [[nodiscard]] bool empty() const { return spfOdd.empty(); }

    /// Returns the effective size of this SPF sieve, which is 1 more than the max valid input to this sieve.
    [[nodiscard]] size_t size() const { return spfOdd.size() * 2 + 1; }

    /// Accessor: returns the smallest prime factor for any x (1 <= x <= n).
    [[nodiscard]] T operator[](T x) const { return get(x); }
};

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

/// Generates a sieve of the mobius function, given a SPF sieve.
template <typename SPFSieve> constexpr std::vector<int8_t> mobiusSieve(size_t limit, const SPFSieve &spfs)
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
template <std::integral T> constexpr std::vector<int8_t> mobiusSieve(T limit) { return mobiusSieve(limit, SPF{limit}); }

/// Generates a sieve of squarefree numbers up to a given limit.
constexpr std::vector<bool> squarefreeSieve(size_t limit)
{
    std::vector sieve(limit + 1, true);
    sieve[0] = false;
    size_t const s = isqrt(limit);
    for (size_t i = 2; i <= s; ++i)
        if (sieve[i])
            for (size_t j = i * i; j <= limit; j += i * i)
                sieve[j] = false;
    return sieve;
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

// Function to find smallest primitive root of p
template <integral2 T, typename SPFSieve = std::vector<int>>
constexpr std::common_type_t<int64_t, T> primitiveRoot(const T &p, SPFSieve &&spfs = {})
{
    using Tp = std::common_type_t<int64_t, T>;
    if (p == 2)
        return Tp(1);
    auto phi = p - 1;
    auto pf = factor(phi, std::forward<SPFSieve>(spfs));
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

/// Calculates `Σ k ≥ 0, f(⌊n / (a * k + b)⌋)`.
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

/// Calculates `Σ k ≥ 0, f(⌊n / (a * k + b)⌋)`.
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
} // namespace euler
