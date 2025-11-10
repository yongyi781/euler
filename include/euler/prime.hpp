#pragma once

/// Constexpr replacement of boost's is_small_prime function.
#include "algorithm.hpp"
#include "literals.hpp"
#include "modular_arithmetic.hpp"
#include <algorithm>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <execution>
#include <primesieve.hpp>
#include <ranges>

namespace euler
{
namespace detail
{
inline constexpr std::array<u32, 17> primesTo59{2U,  3U,  5U,  7U,  11U, 13U, 17U, 19U, 23U,
                                                29U, 31U, 37U, 41U, 43U, 47U, 53U, 59U};
inline constexpr uint32_t firstNontrivialPrime = 61;
inline constexpr std::array primeWheel{10U, 2U, 4U, 2U, 4U, 6U, 2U, 6U, 4U, 2U, 4U, 6U, 6U, 2U, 6U,  4U,
                                       2U,  6U, 4U, 6U, 8U, 4U, 2U, 4U, 2U, 4U, 8U, 6U, 4U, 6U, 2U,  4U,
                                       6U,  2U, 6U, 6U, 4U, 2U, 4U, 6U, 2U, 6U, 4U, 2U, 4U, 2U, 10U, 2U};
inline constexpr uint32_t primeWheelPeriod = 210;
inline constexpr size_t primeWheelSize = std::size(primeWheel);

/// Assumption: n <= 227.
constexpr bool isSmallPrime(size_t n)
{
    std::initializer_list<u32> const primesTo227{
        3U,   5U,   7U,   11U,  13U,  17U,  19U,  23U,  29U,  31U,  37U,  41U,  43U,  47U,  53U,  59U,
        61U,  67U,  71U,  73U,  79U,  83U,  89U,  97U,  101U, 103U, 107U, 109U, 113U, 127U, 131U, 137U,
        139U, 149U, 151U, 157U, 163U, 167U, 173U, 179U, 181U, 191U, 193U, 197U, 199U, 211U, 223U, 227U};
    return std::ranges::contains(primesTo227, n);
}

template <integral2 T> constexpr bool checkSmallFactors(const T &n)
{
    constexpr u32 pp1 = 223092870U;
    auto const m1 = (u32)(n % pp1);
    for (u32 const p : {3U, 5U, 7U, 11U, 13U, 17U, 19U, 23U})
        if (m1 % p == 0)
            return false;

    constexpr u32 pp2 = 2756205443U;
    auto const m2 = (u32)(n % pp2);
    for (u32 const p : {29U, 31U, 37U, 41U, 43U, 47U})
        if (m2 % p == 0)
            return false;

    constexpr u32 pp3 = 907383479U;
    auto const m3 = (u32)(n % pp3);
    for (u32 const p : {53U, 59U, 61U, 67U, 71U})
        if (m3 % p == 0)
            return false;

    constexpr u32 pp4 = 4132280413U;
    auto const m4 = (u32)(n % pp4);
    for (u32 const p : {73U, 79U, 83U, 89U, 97U})
        if (m4 % p == 0)
            return false;

    std::array<std::initializer_list<u32>, 6> smallFactors5{{{101U, 103U, 107U, 109U},
                                                             {113U, 127U, 131U, 137U},
                                                             {139U, 149U, 151U, 157U},
                                                             {163U, 167U, 173U, 179U},
                                                             {181U, 191U, 193U, 197U},
                                                             {199U, 211U, 223U, 227U}}};
    constexpr std::array pp5{121330189U,
                             113U * 127U * 131U * 137U,
                             139U * 149U * 151U * 157U,
                             163U * 167U * 173U * 179U,
                             181U * 191U * 193U * 197U,
                             199U * 211U * 223U * 227U};

    for (std::size_t k = 0; k < pp5.size(); ++k)
    {
        auto const m5 = (u32)(n % pp5[k]);
        for (auto p : smallFactors5[k])
            if (m5 % p == 0)
                return false;
    }
    return true;
}

/// Deterministic Miller test.
template <integral2 T> constexpr bool miller(const T &n, std::initializer_list<u64> witnesses)
{
    T const m = n - 1;
    size_t const k = boost::multiprecision::lsb(m);
    T const q = m >> k;

    T y;
    for (u64 x : witnesses)
    {
        if (x >= n)
        {
            x %= (u64)n;
            if (x == 0)
                return true;
        }
        y = powmSafe(T(x), q, n);
        size_t j = 0;
        while (true)
        {
            if (y == m)
                break;
            if (y == 1)
            {
                if (j == 0)
                    break;
                return false;
            }
            if (++j == k)
                return false;
            y = modmul(y, y, n);
        }
    }
    return true;
}

/// Probabilistic Miller-Rabin test.
template <integral2 T> bool millerRabin(const T &n, size_t trials)
{
    T const m = n - 1;
    size_t const k = boost::multiprecision::lsb(m);
    T const q = m >> k;

    static std::mt19937_64 gen;
    boost::multiprecision::uniform_int_distribution<T> const dist(2, n - 2);

    T x, y;
    for (std::size_t i = 0; i < trials; ++i)
    {
        x = dist(gen);
        y = powmSafe(x, q, n);
        size_t j = 0;
        while (true)
        {
            if (y == m)
                break;
            if (y == 1)
            {
                if (j == 0)
                    break;
                return false;
            }
            if (++j == k)
                return false;
            y = modmul(y, y, n);
        }
    }
    return true;
}
} // namespace detail

/// Returns whether the given number is prime.
/// @param trials The number of Miller-Rabin trials to perform. Only used if n is at least
///               3825123056546413051.
template <integral2 T> constexpr bool isPrime(const T &n, size_t trials = 8)
{
    if (n < 2)
        return false;
    if (boost::multiprecision::bit_test(n, 0) == 0)
        return n == 2; // n is even
    if (n <= 227)
        return euler::detail::isSmallPrime((size_t)n);
    if (!euler::detail::checkSmallFactors(n))
        return false;
    if (n < 52441)
        return true;
    if (n < 341531)
        return euler::detail::miller(n, {9345883071009581737_u64});
    if (n < 1050535501)
        return euler::detail::miller(n, {336781006125_u64, 9639812373923155_u64});
    if (n < 350269456337)
        return euler::detail::miller(n, {4230279247111683200_u64, 14694767155120705706_u64, 16641139526367750375_u64});
    if (n < 55245642489451)
        return euler::detail::miller(n,
                                     {2_u64, 141889084524735_u64, 1199124725622454117_u64, 11096072698276303650_u64});
    if (n < 7999252175582851)
        return euler::detail::miller(
            n, {2_u64, 4130806001517_u64, 149795463772692060_u64, 186635894390467037_u64, 3967304179347715805_u64});
    if (n < 585226005592931977)
        return euler::detail::miller(n, {2_u64, 123635709730000_u64, 9233062284813009_u64, 43835965440333360_u64,
                                         761179012939631437_u64, 1263739024124850375_u64});
    if constexpr (sizeof(T) <= 8)
        return euler::detail::miller(n, {2_u64, 325_u64, 9375_u64, 28178_u64, 450775_u64, 9780504_u64, 1795265022_u64});
    else
        return euler::detail::millerRabin(n, trials);
}

/// Returns whether the given number is prime, using GMP.
template <> inline bool isPrime<mpz_int>(const mpz_int &n, size_t trials)
{
    return mpz_probab_prime_p(n.backend().data(), trials);
}

/// Removes all factors of p from a number, and returns the number of factors removed. The `knownDivides` template
/// parameter is for performance optimization, skipping the first modulus check if it is known `p` divides `n` in
/// advance.
template <bool KnownDivides = false, typename T, typename U> constexpr int removeFactors(T &n, const U &p)
{
    assert(n != 0);
    assert(p >= 2);

    /// Specialization for GMP integers.
    if constexpr (std::same_as<T, mpz_int> && std::same_as<U, mpz_int>)
        return mpz_remove((mpz_ptr)n.backend().data(), (mpz_srcptr)n.backend().data(), (mpz_srcptr)p.backend().data());
    if constexpr (std::integral<T>)
    {
        if (p == 2)
        {
            int e = std::countr_zero(std::make_unsigned_t<T>(n));
            n >>= e;
            return e;
        }
    }
    if constexpr (KnownDivides)
    {
        int result = 1;
        n /= p;
        while (n % p == 0)
        {
            ++result;
            n /= p;
        }
        return result;
    }
    else
    {
        int result = 0;
        while (n % p == 0)
        {
            ++result;
            n /= p;
        }
        return result;
    }
}

/// Calculates the valuation of n with respect to a prime p.
template <bool KnownDivides = false, typename T, typename U>
constexpr int valuation(T n, const U &p) // Pass by value intentional
{
    return removeFactors<KnownDivides>(n, p);
}

/// Removes all factors of p from a number, and returns the result. The `knownDivides` template
/// parameter is for performance optimization, skipping the first modulus check if it is known `p` divides `n` in
/// advance.
template <bool KnownDivides = false, typename T, typename U> constexpr T removedFactors(T n, const U &p)
{
    removeFactors<KnownDivides>(n, p);
    return n;
}

/// Calculates the p-adic valuation of n!.
template <integral2 T, integral2 U> constexpr T factorialValuation(T n, const U &p) // Pass by value intentional
{
    T res = 0;
    while (n > 1)
        res += (n /= p);
    return res;
}

/// Calculates the p-adic valuation of binomial(n,k).
template <integral2 T, integral2 U, integral2 V>
constexpr std::common_type_t<T, U> binomValuation(const T &n, const U &k, const V &p)
{
    return factorialValuation(n, p) - factorialValuation(n - k, p) - factorialValuation(k, p);
}

template <integral2 T, bool PrimalityCheck = !std::same_as<T, mpz_int>>
constexpr T smallestPrimeFactor(const T &num, const T &start = 2)
{
    if (start < T(detail::firstNontrivialPrime))
        for (auto p : detail::primesTo59)
            if (num % p == 0)
                return p;
    if constexpr (PrimalityCheck)
        if (num > 100'000'000 && isPrime(num))
            return num;
    size_t c = 14;
    for (T p = start < 61 ? 61 : 61 + T((start - 61) / detail::primeWheelPeriod * detail::primeWheelPeriod);
         p * p <= num; p += detail::primeWheel[c], c = (c + 1 == detail::primeWheelSize ? 0 : c + 1))
        if (num % p == 0)
            return p;
    return num;
}

/// Calculates Euler's totient function, given a factorization.
template <integral2 T, std::ranges::range Range>
constexpr T totient(T n, Range &&factorization) // Pass by value intentional
{
    for (auto &&[p, e] : std::forward<Range>(factorization))
        n -= n / p;
    return n;
}

/// @brief Generates primes within a range.
/// @param start Inclusive lower bound.
/// @param stop Inclusive upper bound. If this is not specified, then use start as stop.
/// @return Primes from start to stop inclusive.
template <typename T = int64_t> constexpr std::vector<T> primeRange(int64_t start, int64_t stop)
{
    std::vector<T> result;
    if (std::is_constant_evaluated() && stop - start < 7000)
    {
        for (int64_t i = start; i <= stop; ++i)
            if (isPrime(i))
                result.push_back(i);
    }
    else
    {
        primesieve::generate_primes(start, stop, &result);
    }
    return result;
}

/// @brief Generates primes within a range, starting at 2.
/// @param stop Inclusive upper bound. If this is not specified, then use start as stop.
/// @return Primes from start to stop inclusive.
template <typename T = int64_t> constexpr std::vector<T> primeRange(int64_t stop) { return primeRange<T>(2, stop); }

/// Sieves numbers of the form n² + 1. Sends triples (n, p, e) to f, over all primes p ^ e
/// exactly dividing n² + 1 for n between 1 and limit. For some reason, calling this function with
/// an unsigned argument gives a performance boost.
template <execution_policy Exec, integral2 T, std::invocable<T, T, int> Fun>
void sieveSquarePlusOne(Exec &&exec, const T &limit, Fun f)
{
    static constexpr auto callback = [](T x, T &y, T p, Fun f, bool known = false) {
        int e = 0;
        if (known)
        {
            y /= p;
            ++e;
        }
        e += removeFactors(y, p);
        if (e > 0)
            f(x, p, e);
    };
    auto s = mapv(std::forward<Exec>(exec), range(T(0), limit), [&](T n) {
        T y = n * n + 1;
        callback(n, y, 2, f);
        callback(n, y, 5, f);
        callback(n, y, 13, f);
        callback(n, y, 17, f);
        callback(n, y, 29, f);
        callback(n, y, 37, f);
        callback(n, y, 41, f);
        callback(n, y, 53, f);
        callback(n, y, 61, f);
        callback(n, y, 73, f);
        callback(n, y, 89, f);
        callback(n, y, 97, f);
        callback(n, y, 101, f);
        callback(n, y, 109, f);
        callback(n, y, 113, f);
        callback(n, y, 137, f);
        callback(n, y, 149, f);
        callback(n, y, 157, f);
        callback(n, y, 173, f);
        callback(n, y, 181, f);
        callback(n, y, 193, f);
        callback(n, y, 197, f);
        callback(n, y, 229, f);
        callback(n, y, 233, f);
        callback(n, y, 241, f);
        callback(n, y, 257, f);
        callback(n, y, 269, f);
        callback(n, y, 277, f);
        callback(n, y, 281, f);
        callback(n, y, 293, f);
        return y;
    });
    for (T n = 20; n <= limit; ++n)
    {
        auto p = s[n];
        if (p == 1)
            continue;
        // Divide everything that's +/-n mod p by p
        for (T k = n; k <= limit; k += p)
            callback(k, s[k], p, f, true);
        for (T k = p - n; k <= limit; k += p)
            callback(k, s[k], p, f, true);
    }
}

/// Sieves numbers of the form n² + 1. Sends triples (n, p, e) to f, over all primes p ^ e
/// exactly dividing n² + 1 for n between 1 and limit.
template <integral2 T, std::invocable<T, T, int> Fun> void sieveSquarePlusOne(const T &limit, Fun f)
{
    sieveSquarePlusOne(std::execution::seq, limit, f);
}

/// Returns the sum of a function `f` over the integers coprime to the given prime list in the range [1, limit]. The
/// function `f` is passed in as its summatory function `F`. For example, to count the coprimes, use `F = identity`.
template <typename SummatoryFun, typename It, typename T>
constexpr std::remove_cvref_t<std::invoke_result_t<SummatoryFun, T>> sumCoprime(SummatoryFun F, It primeBegin,
                                                                                It primeEnd, T limit)
{
    if (primeBegin == primeEnd)
        return F(limit);
    if (primeBegin + 1 == primeEnd)
        return F(limit) - F(limit / *(primeEnd - 1));
    auto const p = *(primeEnd - 1);
    auto res = sumCoprime(F, primeBegin, std::prev(primeEnd), limit);
    if (limit >= p)
        res -= sumCoprime(std::move(F), primeBegin, std::prev(primeEnd), limit / p);
    return res;
}

/// Returns the number of integers coprime to the given prime list in the range [1, limit].
template <typename It, typename T> constexpr T countCoprime(It primeBegin, It primeEnd, T limit)
{
    return sumCoprime(std::identity{}, primeBegin, primeEnd, limit);
}
} // namespace euler
