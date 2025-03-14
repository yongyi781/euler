#pragma once

/// Constexpr replacement of boost's is_small_prime function.
#include "algorithm.hpp"
#include "modular_arithmetic.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <execution>
#include <primesieve.hpp>
#include <ranges>

inline namespace euler
{
namespace detail
{
constexpr std::array<unsigned, 17> primesTo59{2U,  3U,  5U,  7U,  11U, 13U, 17U, 19U, 23U,
                                              29U, 31U, 37U, 41U, 43U, 47U, 53U, 59U};
constexpr uint32_t firstNontrivialPrime = 61;
constexpr std::array primeWheel{10U, 2U, 4U, 2U, 4U, 6U, 2U, 6U, 4U, 2U, 4U, 6U, 6U, 2U, 6U,  4U,
                                2U,  6U, 4U, 6U, 8U, 4U, 2U, 4U, 2U, 4U, 8U, 6U, 4U, 6U, 2U,  4U,
                                6U,  2U, 6U, 6U, 4U, 2U, 4U, 6U, 2U, 6U, 4U, 2U, 4U, 2U, 10U, 2U};
constexpr uint32_t primeWheelPeriod = 210;
constexpr size_t primeWheelSize = std::size(primeWheel);

/// Assumption: n <= 227.
constexpr bool isSmallPrime(size_t n)
{
    std::initializer_list<unsigned> const primesTo227{
        3U,   5U,   7U,   11U,  13U,  17U,  19U,  23U,  29U,  31U,  37U,  41U,  43U,  47U,  53U,  59U,
        61U,  67U,  71U,  73U,  79U,  83U,  89U,  97U,  101U, 103U, 107U, 109U, 113U, 127U, 131U, 137U,
        139U, 149U, 151U, 157U, 163U, 167U, 173U, 179U, 181U, 191U, 193U, 197U, 199U, 211U, 223U, 227U};
    return std::find(primesTo227.begin(), primesTo227.end(), n) != primesTo227.end();
}

template <integral2 T> constexpr bool checkSmallFactors(const T &n)
{
    constexpr unsigned pp1 = 223092870U;

    auto m1 = (unsigned)(n % pp1);
    for (unsigned const p : {3U, 5U, 7U, 11U, 13U, 17U, 19U, 23U})
        if (m1 % p == 0)
            return false;

    constexpr unsigned pp2 = 2756205443U;

    auto m2 = (unsigned)(n % pp2);
    for (unsigned const p : {29U, 31U, 37U, 41U, 43U, 47U})
        if (m2 % p == 0)
            return false;

    constexpr unsigned pp3 = 907383479U;

    auto m3 = (unsigned)(n % pp3);
    for (unsigned const p : {53U, 59U, 61U, 67U, 71U})
        if (m3 % p == 0)
            return false;

    constexpr unsigned pp4 = 4132280413U;

    auto m4 = (unsigned)(n % pp4);
    for (unsigned const p : {73U, 79U, 83U, 89U, 97U})
        if (m4 % p == 0)
            return false;

    std::array<std::initializer_list<unsigned>, 6> smallFactors5{{{101U, 103U, 107U, 109U},
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
        auto m5 = (unsigned)(n % pp5[k]);
        for (auto p : smallFactors5[k])
            if (m5 % p == 0)
                return false;
    }
    return true;
}

template <integral2 T> constexpr bool millerRabinDeterministic(const T &n, std::initializer_list<int> witnesses)
{
    if (n == 2)
        return true; // Trivial special case.
    if (boost::multiprecision::bit_test(n, 0) == 0)
        return false; // n is even
    if (n <= 227)
        return detail::isSmallPrime((size_t)n);

    if (!detail::checkSmallFactors(n))
        return false;

    T nm1 = n - 1;
    T q = n - 1;
    std::size_t k = boost::multiprecision::lsb(q);
    q >>= k;

    // Execute the trials:
    for (auto &&x : witnesses)
    {
        T y = powmSafe(T(x), q, n);
        size_t j = 0;
        while (true)
        {
            if (y == nm1)
                break;
            if (y == 1)
            {
                if (j == 0)
                    break;
                return false;
            }
            if (++j == k)
                return false;
            y = powmSafe(y, 2, n);
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
    if (n < 1373653)
        return detail::millerRabinDeterministic(n, {2, 3});
    if (n < 9080191)
        return detail::millerRabinDeterministic(n, {31, 73});
    if (n < 4759123141)
        return detail::millerRabinDeterministic(n, {2, 7, 61});
    if (n < 1122004669633)
        return detail::millerRabinDeterministic(n, {2, 13, 23, 1662803});
    if (n < 2152302898747)
        return detail::millerRabinDeterministic(n, {2, 3, 5, 7, 11});
    if (n < 3474749660383)
        return detail::millerRabinDeterministic(n, {2, 3, 5, 7, 11, 13});
    if (n < 341550071728321)
        return detail::millerRabinDeterministic(n, {2, 3, 5, 7, 11, 13, 17});
    if (n < 3825123056546413051)
        return detail::millerRabinDeterministic(n, {2, 3, 5, 7, 11, 13, 17, 19, 23});
    return boost::multiprecision::miller_rabin_test(n, trials);
}

/// Removes all factors of p from a number, and returns the number of factors removed. The `knownDivides` template
/// parameter is for performance optimization, skipping the first modulus check if it is known `p` divides `n` in
/// advance.
template <bool KnownDivides = false, typename Tn, typename Tp> constexpr int removeFactors(Tn &n, const Tp &p)
{
    assert(n != 0);
    assert(p >= 2);

    /// Specialization for GMP integers.
    if constexpr (std::same_as<Tn, mpz_int> && std::same_as<Tp, mpz_int>)
        return mpz_remove((mpz_ptr)n.backend().data(), (mpz_srcptr)n.backend().data(), (mpz_srcptr)p.backend().data());
    if constexpr (std::integral<Tn>)
    {
        if (p == 2)
        {
            int e = std::countr_zero(std::make_unsigned_t<Tn>(n));
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
template <typename Tn, typename Tp> constexpr int valuation(Tn n, const Tp &p) // Pass by value intentional
{
    return removeFactors(n, p);
}

/// @brief Calculates the p-adic valuation of n!.
/// @param n A number.
/// @param p A prime number.
/// @return The p-adic valuation of n!.
template <integral2 Tn, integral2 Tp> constexpr Tn factorialValuation(Tn n, const Tp &p) // Pass by value intentional
{
    Tn total = 0;
    while (n > 1)
        total += (n /= p);
    return total;
}

/// Calculates the p-adic valuation of binomial(n,k).
template <integral2 Tn, integral2 Tk, integral2 Tp>
constexpr std::common_type_t<Tn, Tk> binomValuation(const Tn &n, const Tk &k, const Tp &p)
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

/// @brief Sieves primes up to a limit using the sieve of Eratosthenes.
/// @param limit Inclusive upper bound.
/// @return The sieve.
constexpr std::vector<bool> primeSieve(size_t limit)
{
    std::vector<bool> sieve(limit + 1, false);
    auto sequentialSievePass = [&](int i) {
        if (sieve[i])
        {
            auto jLimit = limit / i;
            for (size_t j = i; j <= jLimit; ++j)
                sieve[i * j] = false;
        }
    };
    for (auto p : detail::primesTo59)
    {
        if (limit < p)
            return sieve;
        sieve[p] = true;
    }
    for (int64_t i = 11, c = 1; i <= limit;
         i += detail::primeWheel[c], c = (c + 1 == detail::primeWheelSize ? 0 : c + 1))
        sieve[i] = true;
    int const sqrtLimit = (int)isqrt(limit);
    for (int i = 11; i <= 31; i += 2)
        sequentialSievePass(i);
    // After i = 31, safe to parallelize.
    if (limit >= 10'000'000)
    {
        for (int i = 37; i <= sqrtLimit; i += 2)
            if (sieve[i])
            {
                int64_t const jLimit = limit / i;
                std::for_each(std::execution::par, counting_iterator((int64_t)i), counting_iterator(jLimit + 1),
                              [&](int64_t j) { sieve[i * j] = false; });
            }
    }
    else
    {
        for (int i = 37; i <= sqrtLimit; i += 2)
            sequentialSievePass(i);
    }
    return sieve;
}

/// @brief Generates primes within a range.
/// @param start Inclusive lower bound.
/// @param stop Inclusive upper bound. If this is not specified, then use start as stop.
/// @return Primes from start to stop inclusive.
template <integral2 T> constexpr std::vector<T> primeRange(T start, T stop)
{
    std::vector<T> result;
    if (std::is_constant_evaluated() && stop - start < 7000)
    {
        for (T i = start; i <= stop; ++i)
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
template <integral2 T> constexpr std::vector<T> primeRange(T stop) { return primeRange(T(2), stop); }

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
