#pragma once

#include "it/base.hpp"
#include "prime.hpp"

namespace euler
{
/// Runs a wavefront algorithm in parallel. Limit is inclusive.
template <std::integral T, typename Fun> void par_wavefront(T begin, T end, Fun f, T r = 2)
{
    assert(begin > 0 && "begin must be > 0");
    assert(r > 1 && "r must be > 1");
    for (; begin <= end; begin *= r)
        tbb::parallel_for(tbb::blocked_range<T>(begin, std::min(end + 1, begin * r)), [&](tbb::blocked_range<T> range) {
            for (T n = range.begin(); n < range.end(); ++n)
                f(n);
        });
}

// Space-optimized structure for smallest prime factors (SPF) up to n.
// It stores SPF only for odd numbers. For any even number, the SPF is 2.
template <std::integral T = int64_t> class SPF
{
    // spf_odd[n] is the smallest prime factor of 2*n + 1.
    std::vector<std::make_unsigned_t<half_integer_t<T>>> spf_odd;
    // inv_odd[n] is the bit inverse of 2*n + 1.
    std::vector<std::make_unsigned_t<T>> inv_odd;

  public:
    // Only need to store integer half the size of the input!
    using unsigned_type = std::make_unsigned_t<T>;
    using half_type = std::make_unsigned_t<half_integer_t<T>>;

    SPF() = default;
    explicit SPF(T n) : spf_odd((n + 1) / 2, 0), inv_odd((isqrt(n) + 1) / 2)
    {
        auto const small_primes = primeRange<half_type>(3, isqrt(n));
        tbb::parallel_for(
            tbb::blocked_range(0UZ, spf_odd.size(), 65536UZ),
            [&](tbb::blocked_range<size_t> r) {
                size_t const seg_end = 2 * r.end();
                for (size_t const p : small_primes)
                {
                    if (p * p >= seg_end)
                        break;
                    size_t const min_j = std::max(p * p / 2, (2 * r.begin() + p) / (2 * p) * p + p / 2);
                    for (size_t j = min_j; j < r.end(); j += p)
                        if (spf_odd[j] == 0)
                            spf_odd[j] = p;
                }
            },
            tbb::simple_partitioner{});
        for (size_t i = 0; i < inv_odd.size(); ++i)
            inv_odd[i] = bitInverse(unsigned_type(2 * i + 1));
    }

    /// Returns whether the SPF sieve is empty.
    [[nodiscard]] constexpr bool empty() const noexcept { return spf_odd.empty(); }

    /// Returns the effective size of this SPF sieve, which is 1 more than the max valid input to this sieve.
    [[nodiscard]] constexpr size_t size() const noexcept { return spf_odd.size() * 2 + 1; }

    /// Returns the smallest prime factor for any x (1 ≤ x ≤ size()).
    [[nodiscard]] T operator[](T n) const
    {
        if (n < 2)
            return 0;
        if (n % 2 == 0)
            return 2;
        return spf_odd[n / 2] == 0 ? n : spf_odd[n / 2];
    }

    /// Returns whether the given number is prime. Requires `1 ≤ n ≤ size()`.
    [[nodiscard]] bool isPrime(T n) const { return (*this)[n] == n; }

    /// Returns multiplicative inverse of n in T. Requirement: n is odd and <= sqrt(N).
    [[nodiscard]] T inv(T n) const { return inv_odd[n / 2]; }

    /// Divides n by p in-place. Requirement: n and p are non-negative, p is odd and <= sqrt(N).
    template <std::integral U>
        requires(sizeof(U) <= sizeof(T))
    void div(U &n, T p) const
    {
        n = unsigned_type(n) * inv(p);
    }

    /// Removes the smallest prime factor `p^e` from a number `n < size()`, and return the `(p, e)` pair.
    template <std::integral U>
        requires(sizeof(U) <= sizeof(T))
    PrimePower<U> remove(U &n) const
    {
        if (n % 2 == 0)
        {
            int const e = std::countr_zero(std::make_unsigned_t<U>(n));
            n >>= e;
            return PrimePower<U>{2, e};
        }
        PrimePower<U> res{(*this)[n], 1};
        if (n == res.first)
        {
            n = 1;
            return res;
        }
        for (div(n, res.first); (*this)[n] == res.first; div(n, res.first), ++res.second)
            ;
        return res;
    }

    /// Removes the smallest prime factor `p^e` from a number `n < size()`, and return `p^e`.
    template <std::integral U>
        requires(sizeof(U) <= sizeof(T))
    U extract(U &n) const
    {
        if (n % 2 == 0)
        {
            int const e = std::countr_zero(std::make_unsigned_t<U>(n));
            n >>= e;
            return U(1) << e;
        }
        T const p = (*this)[n];
        if (n == p)
        {
            n = 1;
            return p;
        }
        U res = p;
        for (div(n, p); (*this)[n] == p; div(n, p), res *= p)
            ;
        return res;
    }

    /// Returns the largest prime factor of `n`.
    [[nodiscard]] T lpf(T n) const
    {
        T res = 0;
        while (n != 1)
            res = remove(n).first;
        return res;
    }

    /// Returns the radical of the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] T radical(T n) const
    {
        T res = 1;
        while (n != 1)
            res *= remove(n).first;
        return res;
    }

    /// Returns the mobius function for the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] int mobius(T n) const
    {
        int res = 1;
        if (n % 2 == 0)
        {
            if (n % 4 == 0)
                return 0;
            n /= 2;
            res = -res;
        }
        while (n != 1)
        {
            auto const p = (*this)[n];
            if (n == p)
                return -res;
            div(n, p);
            if ((*this)[n] == p)
                return 0;
            res = -res;
        }
        return res;
    }

    /// Returns whether `n` is squarefree. Requires 1 ≤ n ≤ size().
    [[nodiscard]] bool sqfree(T n) const
    {
        if (n % 2 == 0)
        {
            if (n % 4 == 0)
                return false;
            n /= 2;
        }
        while (n != 1)
        {
            auto const p = (*this)[n];
            if (n == p)
                return true;
            div(n, p);
            if ((*this)[n] == p)
                return false;
        }
        return true;
    }

    [[nodiscard]] T totient(T n) const
    {
        T res = n;
        while (n != 1)
            res -= (double)res / remove(n).first; // This works for all reasonably sized N.
        return res;
    }

    [[nodiscard]] T countDivisors(T n) const
    {
        T res = 1;
        while (n != 1)
            res *= remove(n).second + 1;
        return res;
    }

    /// Returns the number of distinct prime factors of the given number. Requires 1 ≤ n ≤ size().
    [[nodiscard]] u32 omega(T n) const
    {
        u32 res = 0;
        while (n != 1)
        {
            ++res;
            remove(n);
        }
        return res;
    }

    /// Returns the number of prime factors of the given number, with multiplicity. Requires 1 ≤ n ≤ size().
    [[nodiscard]] u32 Omega(T n) const
    {
        u32 res = 0;
        while (n != 1)
            res += remove(n).second;
        return res;
    }

    /// Returns the number of integers coprime to the given prime list in the range [1, limit].
    template <typename Tk> [[nodiscard]] T countCoprime(Tk k, T limit) const
    {
        return sumCoprime([](auto &&) { return T(1); }, std::identity{}, k, limit);
    }

    /// Decomposes a number into its squarefree part and its square part.
    ///
    /// @param n The number to decompose.
    /// @return A pair of (s, t) such that n = s^2 * t with t squarefree.
    [[nodiscard]] std::pair<T, T> sqfreeDecompose(T n) const
    {
        std::pair<T, T> res{1, 1};
        if (n % 2 == 0)
        {
            int const e = std::countr_zero(unsigned_type(n));
            n >>= e;
            res.first <<= e / 2;
            res.second <<= e % 2;
        }
        while (n != 1)
        {
            T const p = (*this)[n];
            if (n == p)
            {
                n = 1;
                res.second *= p;
                break;
            }
            div(n, p);
            if ((*this)[n] == p)
            {
                res.first *= p;
                div(n, p);
            }
            else
            {
                res.second *= p;
            }
        }
        return res;
    }

    // ==== Sieves ====

    /// Creates a sieve from a multiplicative function in O(n) time.
    template <typename Fun> [[nodiscard]] auto sieve(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T, int>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1);
        res[1] = 1;
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
            {
                int const e = std::countr_zero(unsigned_type(n));
                res[n] = res[n >> e] * f(T(2), e);
            }
            else
            {
                T x = n;
                auto const [p, e] = remove(x);
                res[n] = res[x] * f(p, e);
            }
        });
        return res;
    }

    /// Creates a sieve from an additive function in O(n log log n) time.
    template <typename Fun> [[nodiscard]] auto sieveAdditive(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T, int>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1);
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
            {
                int const e = std::countr_zero(unsigned_type(n));
                res[n] = res[n >> e] + f(T(2), e);
            }
            else
            {
                T x = n;
                auto const [p, e] = remove(x);
                res[n] = res[x] + f(p, e);
            }
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
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
                res[n] = res[n >> 1] * f(T(2));
            else
            {
                T const p = (*this)[n];
                res[n] = res[n / p] * f(p);
            }
        });
        return res;
    }

    /// Creates a sieve from a completely additive function in O(n) time.
    template <typename Fun> [[nodiscard]] auto sieveCompletelyAdditive(Fun f, T limit) const
    {
        using Tp = std::remove_cvref_t<std::invoke_result_t<Fun, T>>;
        assert(limit < size());
        std::vector<Tp> res(limit + 1);
        res[1] = 0;
        for (T start = 2; start <= limit; start <<= 1)
            tbb::parallel_for(tbb::blocked_range<T>(start, std::min(limit + 1, start << 1)),
                              [&](tbb::blocked_range<T> r) {
                                  for (T n = r.begin(); n < r.end(); ++n)
                                  {
                                      if (n % 2 == 0)
                                          res[n] = res[n >> 1] + f(T(2));
                                      else
                                      {
                                          T const p = (*this)[n];
                                          res[n] = res[n / p] + f(p);
                                      }
                                  }
                              });
        return res;
    }

    /// Sieve for largest prime factors in O(n log log n) time.
    [[nodiscard]] std::vector<T> lpfSieve(T limit) const
    {
        assert(limit < size());
        std::vector<T> res(limit + 1);
        res[2] = 2;
        res[3] = 3;
        for (T start = 4; start <= limit; start <<= 1)
            tbb::parallel_for(tbb::blocked_range<T>(start, std::min(limit + 1, start << 1)),
                              [&](tbb::blocked_range<T> r) {
                                  for (T n = r.begin(); n < r.end(); ++n)
                                  {
                                      if (n % 2 == 0)
                                          res[n] = res[n >> 1];
                                      else
                                      {
                                          T const p = (*this)[n];
                                          if (p == n)
                                              res[n] = p;
                                          else
                                          {
                                              T x = n;
                                              div(x, p);
                                              res[n] = max(p, res[x]);
                                          }
                                      }
                                  }
                              });
        return res;
    }

    /// Sieve for Euler's totient function in O(n) time.
    [[nodiscard]] std::vector<T> totientSieve(T limit) const
    {
        assert(limit < size());
        std::vector<T> res(limit + 1);
        res[1] = 1;
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
                res[n] = res[n >> 1] << ((~n & 2) >> 1);
            else
            {
                T const p = (*this)[n];
                if (p == n)
                    res[n] = n - 1;
                else
                {
                    T x = n;
                    div(x, p);
                    res[n] = res[x] * (p - (x % p != 0));
                }
            }
        });
        return res;
    }

    /// Sieve for the divisor counting function in O(n) time.
    template <typename U = half_integer_t<T>> [[nodiscard]] std::vector<U> divisorCountSieve(T limit) const
    {
        return sieve([&](T, int e) { return U(e + 1); }, limit);
    }

    /// Sieve for the σ₁ function, the divisor sum function, in O(n) time.
    template <typename U = T> [[nodiscard]] std::vector<U> divisorSumSieve(T limit) const
    {
        assert(limit < size());
        std::vector<U> res(limit + 1);
        res[1] = 1;
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
            {
                int const e = std::countr_zero(unsigned_type(n));
                res[n] = res[n >> e] * ((U(1) << (e + 1)) - 1);
            }
            else
            {
                T const p = (*this)[n];
                U s = 1 + p;
                if (p == n)
                    res[n] = s;
                else
                {
                    T x = n;
                    div(x, p);
                    for (U q = p; (*this)[x] == p; div(x, p), s += (q *= p))
                        ;
                    res[n] = res[x] * s;
                }
            }
        });
        return res;
    }

    /// Sieve for the σ₂ function, the divisor square sum function, in O(n) time.
    template <typename U = T> [[nodiscard]] std::vector<U> divisorSum2Sieve(T limit) const
    {
        assert(limit < size());
        std::vector<U> res(limit + 1);
        res[1] = 1;
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
            {
                int const e = std::countr_zero(unsigned_type(n));
                res[n] = res[n >> e] * (((U(1) << 2 * (e + 1)) - 1) / 3);
            }
            else
            {
                T const p = (*this)[n];
                T const pp = p * p;
                U s = 1 + pp;
                if (p == n)
                    res[n] = s;
                else
                {
                    T x = n;
                    div(x, p);
                    for (U q = pp; (*this)[x] == p; div(x, p), s += (q *= pp))
                        ;
                    res[n] = res[x] * s;
                }
            }
        });
        return res;
    }

    /// Sieve for the Möbius function in O(n) time.
    [[nodiscard]] std::vector<i8> mobiusSieve(T limit) const
    {
        assert(limit < size());
        std::vector<i8> res(limit + 1);
        res[1] = 1;
        par_wavefront(T(2), limit, [&](T n) {
            if (n % 2 == 0)
                res[n] = n % 4 == 0 ? 0 : -res[n >> 1];
            else
            {
                T const p = (*this)[n];
                if (p == n)
                    res[n] = -1;
                else
                {
                    T x = n;
                    div(x, p);
                    res[n] = (*this)[x] == p ? 0 : -res[x];
                }
            }
        });
        return res;
    }

    /// Sieve for the Liouville λ function in O(n) time.
    [[nodiscard]] std::vector<i8> liouvilleSieve(T limit) const
    {
        return sieveCompletelyMultiplicative([](T) -> i8 { return -1; }, limit);
    }

    /// Sieve for the ω function, the number of distinct prime factors of a number in O(n) time.
    [[nodiscard]] std::vector<u8> omegaSieve(T limit) const
    {
        return sieveAdditive([](T, int) -> u8 { return 1; }, limit);
    }

    /// Sieve for the Ω function, the number of prime factors of a number in O(n) time.
    [[nodiscard]] std::vector<u8> OmegaSieve(T limit) const
    {
        return sieveCompletelyAdditive([](T) -> u8 { return 1; }, limit);
    }
};

/// Sieve for Euler's totient function in O(n log log n) time.
template <std::integral T> std::vector<T> totientSieve(T limit) { return SPF{limit}.totientSieve(limit); }

/// Sieve for the divisor counting function.
template <std::integral T> std::vector<half_integer_t<T>> divisorCountSieve(T limit)
{
    return SPF{limit}.divisorCountSieve(limit);
}

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
} // namespace euler
