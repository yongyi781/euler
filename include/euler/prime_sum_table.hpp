#pragma once

#include "FenwickTree.hpp"
#include "floors_array.hpp"
#include "it/primes.hpp"
#include "libdivide.h"

namespace euler
{
/// Returns a floors array of values `(1 ≤ p ≤ k, p prime) * f(p)` for `k` of the form `⌊limit / i⌋`. Here, `f` must be
/// a completely multiplicative function and `F` must be the summatory function of `f`.
template <typename Fun, typename SummatoryFun> auto primeSumTable(size_t N, Fun f, SummatoryFun F, double alpha = 0.2)
{
    using T = std::remove_cvref_t<
        std::common_type_t<std::invoke_result_t<Fun, size_t>, std::invoke_result_t<SummatoryFun, size_t>>>;
    size_t const s_N = isqrt(N);
    size_t const fen_sz = std::max(s_N + 1, (size_t)(alpha * std::pow(N / log(N), 2.0 / 3)));
    std::vector<bool> is_prime(fen_sz / 2 + 1, true);
    is_prime[0] = false;
    T const f2 = T(f(2));
    FenwickTree<T> S(fen_sz / 2 + 1, [&](size_t i) { return i == 0 ? T{} : i == 1 ? f2 + f(3) : f(2 * i + 1); });
    constexpr auto idx = [](size_t k) { return (k - 1) >> 1; };
    floors_array<T> res(N);
    tbb::parallel_for(1UZ, N / fen_sz + 1, [&](u32 i) {
        size_t const n = res.quotient(i);
        res.down(i) = F(n) - F(1) - f2 * (F(n / 2) - F(1));
    });
    size_t mid_p = isqrt(fen_sz);
    if (mid_p % 2 == 0)
        --mid_p;

    std::vector<T> Sv(fen_sz / 2 + 1);
    size_t i_frozen = 0;
    T running_prime_sum{};

    auto const freeze_step = [&] {
        if (is_prime[i_frozen])
        {
            size_t const num = 2 * i_frozen + 1;
            running_prime_sum += (num == 1) ? T{} : (num == 3 ? f2 + f(3) : f(num));
        }
        Sv[i_frozen] = running_prime_sum;
    };

    auto const update_step = [&](size_t p, auto &&sum) {
        T const fp = T(f(p));
        T const prev = p == 3 ? f2 : sum(p - 1);
        size_t const pp = p * p;
        size_t const mid_i = (N - 1) / (p * fen_sz);
        size_t const max_i = N / std::max(fen_sz, pp);
        for (size_t i = 1, j = p; i <= mid_i; ++i, j += p)
            res.down(i) -= fp * (res.down(j) - prev);
        libdivide::divider const fast_p(p);
        if (max_i - mid_i > 2048)
            tbb::parallel_for(mid_i + 1, max_i + 1,
                              [&](size_t i) { res.down(i) -= fp * (sum(res.quotient(i) / fast_p) - prev); });
        else
            for (size_t i = mid_i + 1; i <= max_i; ++i)
                res.down(i) -= fp * (sum(res.quotient(i) / fast_p) - prev);
    };

    for (size_t p = 3; p <= mid_p; p += 2)
    {
        if (!is_prime[p / 2])
            continue;
        size_t const max_frozen = std::min(Sv.size(), idx(p * p));
        for (; i_frozen < max_frozen; ++i_frozen)
            freeze_step();
        update_step(p, [&](size_t k) { return k < i_frozen ? Sv[idx(k)] : S.sum(idx(k)); });
        for (size_t k = p * p; k <= fen_sz; k += 2 * p)
            if (is_prime[k / 2])
            {
                is_prime[k / 2] = false;
                S.add(idx(k), -T(f(k)));
            }
    }
    for (; i_frozen < Sv.size(); ++i_frozen)
        freeze_step();
    for (size_t p = mid_p + 2; p <= s_N; p += 2)
    {
        if (!is_prime[p / 2])
            continue;
        update_step(p, [&](size_t k) { return Sv[idx(k)]; });
    }
    res.up(2) = f2;
    tbb::parallel_for(3UZ, res.up().size(), [&](u32 n) { res.up(n) = Sv[idx(n)]; });
    tbb::parallel_for(N / fen_sz + 1, res.down().size(), [&](u32 i) {
        size_t const n = res.quotient(i);
        res.down(i) = Sv[idx(n)];
    });
    return res;
}

/// Returns a floors array of values `(1 ≤ p ≤ k, p prime) * p` for `k` of the form `⌊limit / i⌋`.
template <typename T = u64> constexpr floors_array<T> primeSumTable(size_t N)
{
    return primeSumTable(N, [](size_t n) -> T { return n; }, [](size_t n) -> T { return sumId<T>(n); });
}

/// Returns a floors array of values `#(1 ≤ p ≤ k, p prime)` for `k` of the form `⌊limit / i⌋`.
template <typename T = u64> constexpr floors_array<T> primePiTable(size_t N)
{
    return primeSumTable(N, [](size_t) -> T { return 1; }, [](size_t n) -> T { return n; });
}

/// Calculates `(1 ≤ p ≤ limit, p prime) * f(p)`. Here, `f` must be a completely multiplicative function and `F` must be
/// the summatory function of `f`.
template <std::invocable<size_t> Fun, std::invocable<size_t> SummatoryFun>
constexpr auto primeSum(size_t N, Fun f, SummatoryFun F)
{
    return primeSumTable(N, std::move(f), std::move(F))[N];
}

/// Calculates `(1 ≤ p ≤ limit, p prime) * p`.
template <typename T = u64> constexpr T primeSum(size_t N)
{
    return primeSum(N, [](size_t n) -> T { return n; }, [](auto &&n) -> T { return sumId<T>(n); });
}

/// Returns a list of pairs `(exp, c)` indicating that `c` primes have exponent `exp` in the factorization of `n!`.
/// O(n^(3/4)). Sublinear version of `factorFactorial`.
template <std::integral T> inline std::vector<std::pair<T, T>> factorialExponents(T N)
{
    T const s = isqrt(N);
    std::vector<std::pair<T, T>> res;
    res.reserve(s + N / (s + 1));
    it::primes(2, s)([&](u64 p) { res.emplace_back(factorialValuation(N, p), 1); });
    auto const S = primePiTable<T>(N);
    for (T i = N / (s + 1); i >= 1; --i)
    {
        T const c = S.down(i) - (i + 1 < S.down().size() ? S.down(i + 1) : S.up(N / (i + 1)));
        if (c > 0)
            res.emplace_back(i, c);
    }
    return res;
}
} // namespace euler
