#pragma once

#include "it/base.hpp"
#include "it/factor.hpp"
#include "prime.hpp"

inline namespace euler
{
// Space-optimized structure for smallest prime factors (SPF) up to n.
// It stores SPF only for odd numbers. For any even number, the SPF is 2.
template <std::integral T = int64_t> class SPF
{
  public:
    // Only need to store integer half the size of the input!
    using half_integer_type = std::make_unsigned_t<half_integer_t<T>>;

    SPF() = default;
    explicit SPF(T n) : spfOdd((n + 1) / 2, 0), smallPrimes(primeRange<half_integer_type>(3, isqrt(n)))
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
            for (half_integer_type const p : smallPrimes)
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

    /// Returns whether the SPF sieve is empty.
    [[nodiscard]] bool empty() const noexcept { return spfOdd.empty(); }

    /// Returns the effective size of this SPF sieve, which is 1 more than the max valid input to this sieve.
    [[nodiscard]] size_t size() const noexcept { return spfOdd.size() * 2 + 1; }

    /// Returns the smallest prime factor for any x (1 ≤ x ≤ size()).
    [[nodiscard]] T operator[](T n) const
    {
        if (n < 2)
            return 0;
        if (n % 2 == 0)
            return 2;
        return spfOdd[n / 2] == 0 ? n : spfOdd[n / 2];
    }

    /// Returns whether the given number is prime. Requires 1 ≤ n ≤ size().
    [[nodiscard]] bool isPrime(T n) const { return (*this)[n] == n; }

    /// Returns the mobius function for the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] int mobius(T n) const
    {
        int res = 1;
        if (n % 4 == 0)
            return 0;
        if (n % 2 == 0)
        {
            n /= 2;
            res = -res;
        }
        while (n != 1)
        {
            auto const p = (*this)[n];
            n /= p;
            if ((*this)[n] == p)
                return 0;
            res = -res;
        }
        return res;
    }

    [[nodiscard]] T totient(T n) const
    {
        T res = n;
        if (n % 2 == 0)
        {
            res >>= 1;
            n >>= std::countr_zero(std::make_unsigned_t<T>(n));
        }
        while (n != 1)
        {
            T const p = (*this)[n];
            res -= res / p;
            while ((*this)[n] == p)
                n /= p;
        }
        return res;
    }

    [[nodiscard]] T countDivisors(T n) const
    {
        T res = 1;
        if (n % 2 == 0)
        {
            int const e = std::countr_zero(std::make_unsigned_t<T>(n));
            res *= e + 1;
            n >>= e;
        }
        while (n != 1)
        {
            T const p = (*this)[n];
            int const e = removeFactors<true>(n, p);
            res *= e + 1;
        }
        return res;
    }

    /// Returns the number of distinct prime factors of the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] uint32_t omega(T n) const
    {
        uint32_t res = 0;
        if (n % 2 == 0)
        {
            ++res;
            n >>= std::countr_zero(std::make_unsigned_t<T>(n));
        }
        while (n != 1)
        {
            auto const p = (*this)[n];
            n /= p;
            while ((*this)[n] == p)
                n /= p;
            ++res;
        }
        return res;
    }

    /// Returns the number of prime factors of the given number, with multiplicity. Requires 1 ≤ n ≤ size().
    [[nodiscard]] uint32_t Omega(T n) const
    {
        uint32_t res = 0;
        if (n % 2 == 0)
        {
            int const e = std::countr_zero(std::make_unsigned_t<T>(n));
            res += e;
            n >>= e;
        }
        while (n != 1)
        {
            auto const p = (*this)[n];
            n /= p;
            ++res;
            while ((*this)[n] == p)
            {
                ++res;
                n /= p;
            }
        }
        return res;
    }

    /// Returns the sum of a function `f` over the integers coprime to `k` in the range [1, limit]. The function `f` is
    /// passed in as its summatory function `F`. For example, to count the coprimes, use `F = identity`.
    template <typename SummatoryFun, typename Tk>
    [[nodiscard]] auto sumCoprime(SummatoryFun F, Tk k, T limit) const
        -> std::remove_cvref_t<std::invoke_result_t<SummatoryFun, T>>
    {
        thread_local std::vector<Tk> primes;
        primes.clear();
        it::factor(k, *this).map([&](auto &&t) { return t.first; }).appendTo(primes);
        return euler::sumCoprime(std::move(F), primes.begin(), primes.end(), limit);
    }

    /// Returns the number of integers coprime to the given prime list in the range [1, limit].
    template <typename Tk> [[nodiscard]] T countCoprime(Tk k, T limit) const
    {
        return sumCoprime(std::identity{}, k, limit);
    }

    /// Decomposes a number into its squarefree part and its square part.
    ///
    /// @param n The number to decompose.
    /// @return A pair of (u, v) such that n = u * v^2 with u squarefree.
    [[nodiscard]] std::pair<T, T> sqfreeDecompose(T n) const
    {
        std::pair<T, T> res{1, n};
        if (n % 4 == 0)
        {
            int const e = std::countr_zero(std::make_unsigned_t<T>(n));
            res.first <<= e / 2;
            res.second >>= 2 * (e / 2);
            n >>= e;
        }
        while (n != 1)
        {
            auto const p = (*this)[n];
            T const q = pow(p, removeFactors<true>(n, p) / 2);
            res.first *= q;
            res.second /= q * q;
        }
        return res;
    }

    // ==== Sieves ====

    /// Creates a sieve from a multiplicative function in O(n log log n) time.
    template <typename Fun> [[nodiscard]] auto sieve(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T, int>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1, 1);
        res[0] = 0;
        it::range(3, limit, 2)(std::execution::par, [&](T i) {
            T n = i;
            while (n != 1 && res[i] != 0)
            {
                T const p = (*this)[n];
                int const e = removeFactors<true>(n, p);
                res[i] *= f(p, e);
            }
        });
        it::range(2, limit, 2)(std::execution::par, [&](T i) {
            int const e = std::countr_zero(std::make_unsigned_t<T>(i));
            res[i] = res[i >> e] * f(2, e);
        });
        return res;
    }

    /// Creates a sieve from an additive function in O(n log log n) time.
    template <typename Fun> [[nodiscard]] auto sieveAdditive(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T, int>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1);
        it::range(3, limit, 2)(std::execution::par, [&](T i) {
            T n = i;
            while (n != 1)
            {
                T const p = (*this)[n];
                int const e = removeFactors<true>(n, p);
                res[i] += f(p, e);
            }
        });
        it::range(2, limit, 2)(std::execution::par, [&](T i) {
            int const e = std::countr_zero(std::make_unsigned_t<T>(i));
            res[i] = res[i >> e] + f(2, e);
        });
        return res;
    }

    /// Creates a sieve from a completely multiplicative function in O(n) time.
    template <typename Fun> [[nodiscard]] auto sieveCompletelyMultiplicative(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1);
        res[1] = 1;
        for (T i = 3; i <= limit; i += 2)
        {
            T const p = (*this)[i];
            res[i] = res[i / p] * f(p);
        }
        for (T i = 2; i <= limit; i += 2)
            res[i] = res[i / 2] * f(2);
        return res;
    }

    /// Creates a sieve from a completely additive function in O(n) time.
    template <typename Fun> [[nodiscard]] auto sieveCompletelyAdditive(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1);
        for (T i = 3; i <= limit; i += 2)
        {
            T const p = (*this)[i];
            res[i] = res[i / p] + f(p);
        }
        for (T i = 2; i <= limit; i += 2)
            res[i] = res[i / 2] + f(2);
        return res;
    }

    /// Sieve for the divisor counting function in O(n log log n) time.
    template <typename U = T> [[nodiscard]] std::vector<U> divisorCountSieve(T limit) const
    {
        return sieve([&](T, int e) -> U { return e + 1; }, limit);
    }

    /// Sieve for the σ₁ function, the divisor sum function in O(n log log n) time.
    template <typename U = T> [[nodiscard]] std::vector<U> divisorSumSieve(T limit) const
    {
        std::vector<U> sieve(limit + 1, 1);
        sieve[0] = 0;
        it::range(3, limit, 2)(std::execution::par, [&](T i) {
            T n = i;
            while (n != 1)
            {
                T const p = (*this)[n];
                n /= p;
                T q = p, S = 1 + p;
                while ((*this)[n] == p)
                {
                    q *= p;
                    n /= p;
                    S += q;
                }
                sieve[i] *= S;
            }
        });
        it::range(2, limit, 2)(std::execution::par, [&](T i) {
            auto e = std::countr_zero(std::make_unsigned_t<T>(i));
            sieve[i] = ((T(1) << (e + 1)) - 1) * sieve[i >> e];
        });
        return sieve;
    }

    /// Sieve for the Möbius function in O(n log log n) time.
    [[nodiscard]] std::vector<int8_t> mobiusSieve(T limit) const
    {
        assert(limit < size());
        std::vector<int8_t> res(limit + 1, 1);
        res[0] = 0;
        it::range(3, limit, 2)(std::execution::par, [&](T i) {
            T n = i;
            while (n != 1)
            {
                T const p = (*this)[n];
                n /= p;
                if ((*this)[n] == p)
                {
                    res[i] = 0;
                    break;
                }
                res[i] = -res[i];
            }
        });
        it::range(2, limit, 2)(std::execution::par, [&](T i) { res[i] = i % 4 == 0 ? 0 : -res[i / 2]; });
        return res;
    }

    /// Sieve for the Liouville λ function in O(n) time.
    [[nodiscard]] std::vector<int8_t> liouvilleSieve(T limit) const
    {
        return sieveCompletelyMultiplicative([](T) -> int8_t { return -1; }, limit);
    }

    /// Sieve for the ω function, the number of distinct prime factors of a number in O(n log log n) time.
    [[nodiscard]] std::vector<uint8_t> omegaSieve(T limit) const
    {
        return sieveAdditive([](T, int) -> uint8_t { return 1; }, limit);
    }

    /// Sieve for the Ω function, the number of prime factors of a number in O(n) time.
    [[nodiscard]] std::vector<uint8_t> OmegaSieve(T limit) const
    {
        return sieveCompletelyAdditive([](T) -> uint8_t { return 1; }, limit);
    }

  private:
    // spfOdd[i] holds the smallest prime factor for number (2*i + 1).
    // Index 0 corresponds to 1 (unused), index 1 to 3, index 2 to 5, etc.
    std::vector<half_integer_type> spfOdd;
    std::vector<half_integer_type> smallPrimes; // Odd primes up to sqrt(n).
};

/// Sieve for the divisor counting function.
template <std::integral T> std::vector<T> divisorCountSieve(T limit) { return SPF{limit}.divisorCountSieve(limit); }

/// Sieve for the σ₁ function, the divisor sum function.
template <std::integral T> std::vector<T> divisorSumSieve(T limit) { return SPF{limit}.divisorSumSieve(limit); }

/// Sieve for the ω function, the number of distinct prime factors of a number.
template <std::integral T> std::vector<uint8_t> omegaSieve(T limit) { return SPF{limit}.omegaSieve(limit); }

/// Sieve for the Möbius function.
template <std::integral T> std::vector<int8_t> mobiusSieve(T limit) { return SPF{limit}.mobiusSieve(limit); }

/// Sieve for the Liouville λ function.
template <std::integral T> std::vector<int8_t> liouvilleSieve(T limit) { return SPF{limit}.liouvilleSieve(limit); }

/// Sieve for the Ω function.
template <std::integral T> std::vector<uint8_t> OmegaSieve(T limit) { return SPF{limit}.OmegaSieve(limit); }

// ==== Sieves that don't use SPF ====

/// Sieve for the totient function. O(n).
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

/// Generates a sieve of squarefree numbers up to a given limit.
template <typename T = bool> std::vector<T> squarefreeSieve(size_t limit)
{
    std::vector<T> sieve(limit + 1, true);
    sieve[0] = false;
    size_t const s = isqrt(limit);
    for (size_t i = 2; i <= s; ++i)
        if (sieve[i])
            for (size_t j = i * i; j <= limit; j += i * i)
                sieve[j] = false;
    return sieve;
}
} // namespace euler
