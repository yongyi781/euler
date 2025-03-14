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
    explicit SPF(T n)
        : spfOdd((n + 1) / 2, 0), smallPrimes(primeRange(half_integer_type(3), half_integer_type(isqrt(n))))
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
    [[nodiscard]] T operator[](T n) const noexcept
    {
        if (n < 2)
            return 0;
        if (n % 2 == 0)
            return 2;
        return spfOdd[n / 2] == 0 ? n : spfOdd[n / 2];
    }

    /// Returns whether the given number is prime. Requires 1 ≤ n ≤ size().
    [[nodiscard]] bool isPrime(T n) const noexcept { return (*this)[n] == n; }

    /// Returns the mobius function for the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] int mobius(T n) const noexcept
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

    /// Returns the number of distinct prime factors of the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] uint32_t omega(T n) const noexcept
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
    [[nodiscard]] uint32_t Omega(T n) const noexcept
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
    template <typename Tk> constexpr T countCoprime(Tk k, T limit) const
    {
        return sumCoprime(std::identity{}, k, limit);
    }

  private:
    // spfOdd[i] holds the smallest prime factor for number (2*i + 1).
    // Index 0 corresponds to 1 (unused), index 1 to 3, index 2 to 5, etc.
    std::vector<half_integer_type> spfOdd;
    std::vector<half_integer_type> smallPrimes; // Odd primes up to sqrt(n).
};

/// Sieve for divisor counts. This is faster than divisorCountSieve2 for limits over 2 million.
template <typename T = int64_t> constexpr std::vector<T> divisorCountSieve(T limit)
{
    std::vector<T> sieve(limit + 1, 1);
    sieve[0] = 0;
    SPF const spfs{T(limit)};
    it::range(3, limit, 2)(std::execution::par, [&](T i) {
        T n = i;
        while (n != 1)
        {
            auto const p = spfs[n];
            int e = 1;
            n /= p;
            while (spfs[n] == p)
            {
                ++e;
                n /= p;
            }
            sieve[i] *= 1 + e;
        }
    });
    it::range(2, limit, 2)(std::execution::par, [&](T i) {
        auto e = std::countr_zero(std::make_unsigned_t<T>(i));
        sieve[i] = (1 + e) * sieve[i >> e];
    });
    return sieve;
}

/// Sieve for the σ₁ function, the divisor sum function.
template <typename T = int64_t> constexpr std::vector<T> divisorSumSieve(T limit)
{
    std::vector<T> sieve(limit + 1, 1);
    sieve[0] = 0;
    SPF const spfs{T(limit)};
    it::range(3, limit, 2)(std::execution::par, [&](T i) {
        T n = i;
        while (n != 1)
        {
            auto const p = spfs[n];
            n /= p;
            T q = p, S = 1 + p;
            while (spfs[n] == p)
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

/// Sieve for the ω function, the number of distinct prime factors of a number.
constexpr std::vector<uint8_t> omegaSieve(size_t limit)
{
    std::vector<uint8_t> ω(limit + 1);
    SPF const spfs{limit};
    for (size_t i = 2; i <= limit; ++i)
    {
        size_t const j = i / spfs[i];
        ω[i] = j % spfs[i] == 0 ? ω[j] : ω[j] + 1;
    }
    return ω;
}

/// Generates a sieve of the mobius function, given a SPF sieve.
template <typename SPFSieve> constexpr std::vector<int8_t> mobiusSieve(size_t limit, const SPFSieve &spfs)
{
    assert(limit <= spfs.size() - 1);
    std::vector<int8_t> μ(limit + 1, 1);
    μ[0] = 0;
    it::range(3, limit, 2)(std::execution::par, [&](size_t i) {
        size_t n = i;
        while (n != 1)
        {
            auto const p = spfs[n];
            n /= p;
            if (spfs[n] == p)
            {
                μ[i] = 0;
                break;
            }
            μ[i] = -μ[i];
        }
    });
    it::range(2, limit, 2)(std::execution::par, [&](size_t i) { μ[i] = i % 4 == 0 ? 0 : -μ[i / 2]; });
    return μ;
}

/// Sieve for the Mobius function.
template <std::integral T> constexpr std::vector<int8_t> mobiusSieve(T limit) { return mobiusSieve(limit, SPF{limit}); }

/// Generates a sieve of squarefree numbers up to a given limit.
template <typename T = bool> constexpr std::vector<T> squarefreeSieve(size_t limit)
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

/// Generates a sieve of the Liouville function, given an SPF sieve.
template <typename SPFSieve> constexpr std::vector<int8_t> liouvilleSieve(size_t limit, const SPFSieve &spfs)
{
    assert(limit <= spfs.size() - 1);
    std::vector<int8_t> λ(limit + 1);
    λ[1] = 1;
    for (size_t i = 3; i <= limit; i += 2)
        λ[i] = -λ[i / spfs[i]];
    for (size_t i = 2; i <= limit; i += 2)
        λ[i] = -λ[i / 2];
    return λ;
}

/// Sieve for the Liouville function.
template <std::integral T> constexpr std::vector<int8_t> liouvilleSieve(T limit)
{
    return liouvilleSieve(limit, SPF{limit});
}

/// Generates a sieve of the Ω function, given a SPF sieve.
template <typename SPFSieve> constexpr std::vector<uint8_t> bigOmegaSieve(size_t limit, const SPFSieve &spfs)
{
    assert(limit <= spfs.size() - 1);
    std::vector<uint8_t> Ω(limit + 1);
    for (size_t i = 3; i <= limit; i += 2)
        Ω[i] = Ω[i / spfs[i]] + 1;
    for (size_t i = 2; i <= limit; i += 2)
        Ω[i] = Ω[i / 2] + 1;
    return Ω;
}

/// Sieve for the Ω function.
template <std::integral T> constexpr std::vector<uint8_t> bigOmegaSieve(T limit)
{
    return bigOmegaSieve(limit, SPF{limit});
}
} // namespace euler
