#pragma once

#include "../SPF.hpp"
#include "../it/primes.hpp"
#include "../math.hpp"
#include "Dirichlet.hpp"

namespace euler::dirichlet
{
// O(1) Dirichlet series come in two flavors, one is a SpecialDirichlet (recommended) and the other is a Dirichlet.

/// 1, the multiplicative identity. f(n) = [n = 1]. Motive = 0.
template <typename T = u64> constexpr auto unit()
{
    return SpecialDirichlet{[](size_t k) -> T { return k == 1; }, [](size_t) -> T { return 1; }};
}

/// 1, the multiplicative identity. f(n) = [n = 1]. Motive = 0.
template <typename T = u64> Dirichlet<T> unit(size_t n)
{
    return {n, [](size_t) -> T { return 1; }};
}

/// ζ(s). f(n) = 1. Motive = [1].
template <typename T = u64> constexpr auto zeta()
{
    return SpecialDirichlet{[](size_t) -> T { return 1; }, [](size_t k) -> T { return k; }};
}

/// ζ(s). f(n) = 1. Motive = [1].
template <typename T = u64> Dirichlet<T> zeta(size_t n) { return {n, std::identity{}}; }

/// ζ(s - 1). f(n) = n. Motive = [p].
template <typename T = u64> constexpr auto id()
{
    return SpecialDirichlet{[](size_t k) -> T { return k; }, [](size_t k) -> T { return sumId<T>(k); }};
}

/// ζ(s - 1). f(n) = n. Motive = [p].
template <typename T = u64> Dirichlet<T> id(size_t n)
{
    return {n, [](size_t k) -> T { return sumId<T>(k); }};
}

/// ζ(s - 2). f(n) = n^2. If `T = ZMod<M>`, currently requires that neither 2 or 3 divide `M`. Motive = [p^2].
template <typename T = u64> constexpr auto id2()
{
    return SpecialDirichlet{[](size_t k) -> T { return T(k) * T(k); }, [](size_t k) -> T { return sumSquares<T>(k); }};
}

/// ζ(s - 2). f(n) = n^2. If `T = ZMod<M>`, currently requires that neither 2 or 3 divide `M`. Motive = [p^2].
template <typename T = u64> Dirichlet<T> id2(size_t n)
{
    // Fancy stuff to avoid overflow or ZMod divisions.
    return {n, [](size_t k) -> T { return sumSquares<T>(k); }};
}

/// ζ(2s). f(n) = [n is square]. Motive = [1] + [-1].
template <typename T = u64> constexpr auto zeta_2s()
{
    return SpecialDirichlet{[](size_t k) -> T { return isSquare(k); }, [](size_t k) -> T { return isqrt(k); }};
}

/// ζ(2s). f(n) = [n is square]. Motive = [1] + [-1].
template <typename T = u64> Dirichlet<T> zeta_2s(size_t n)
{
    return {n, [](size_t k) -> T { return isqrt(k); }};
}

// TODO: Add SpecialDirichlets for these at some point.

/// ζ(s - 3). f(n) = n^3. Motive = [p^3].
template <typename T = u64> Dirichlet<T> id3(size_t n)
{
    return {n, [](size_t k) -> T {
                T const s = sumId<T>(k);
                return s * s;
            }};
}

/// ζ(as). f(n) = [n is a perfect ath power]. Motive = sum([r] for r in a-th roots of unity)
template <typename T = u64> Dirichlet<T> zeta_multiple(int a, size_t n)
{
    return {n, [&](size_t k) -> T { return inth_root(k, a); }};
}

/// χ_(-3)(s). f(n) = (-3|n). 1 if 1 mod 3, -1 if 2 mod 3, 0 otherwise. L-function of the hexagonal lattice and the
/// Eisenstein integers.
template <typename T = i8> Dirichlet<T> chi3(size_t n)
{
    return {n, [&](size_t k) -> T { return k % 3 == 1; }};
}

/// χ_(-4)(s). f(n) = (-4|n). 1 if 1 mod 4, -1 if 3 mod 4, 0 otherwise. Also known as the Dirichlet beta function.
/// Motive = χ_(-4).
template <typename T = i8> Dirichlet<T> chi4(size_t n)
{
    return {n, [&](size_t k) -> T { return k % 4 == 1 || k % 4 == 2; }};
}

/// χ_5(s). f(n) = (5|n).
template <typename T = i8> Dirichlet<T> chi5(size_t n)
{
    return {n, [&](size_t k) -> T { return std::array{0, 1, 0, -1, 0}[k % 5]; }};
}

/// 1 / ζ(2s). f(n) = [n is square] * μ(√n). O(n^(1/2)). Motive = -[1] - [-1].
template <typename T = int> Dirichlet<T> inv_zeta_2s(size_t n)
{
    auto const mu = mobiusSieve(isqrt(n));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](size_t k) -> T { return mertens[isqrt(k)]; }};
}

/// 1 / ζ(as). f(n) = [n is a perfect rth power] * μ(n^(1/a)). O(n^(1/a)). Motive = -∑([r] for r in ath roots of unity).
template <typename T = u64> Dirichlet<T> inv_zeta_multiple(int a, size_t n)
{
    auto const mu = mobiusSieve(inth_root(n, a));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](size_t k) -> T { return mertens[inth_root(k, a)]; }};
}

/// ζ(as - b). f(n) = [n = k^a] * k^b. Requires `a > 1`. Motive = ∑([r*p^(b/a)] for r in ath roots of unity).
template <typename T = u64> Dirichlet<T> zeta_linear(int a, int b, size_t n)
{
    using std::pow;
    assert(a > 1);
    size_t const s = inth_root(n, a);
    auto sieve = range(0, s, [&](size_t k) { return pow(T(k), b); });
    sieve[0] = 0;
    partialSumInPlace(std::execution::par, sieve);
    return {n, [&](size_t k) -> T { return sieve[inth_root(k, a)]; }};
}

/// 1 / ζ(as - b). f(n) = [n = k^a] * μ(k) * k^b. Requires `a > 1`. Motive = -∑([r*p^(b/a)] for r in ath roots of
/// unity).
template <typename T = i64> Dirichlet<T> inv_zeta_linear(int a, int b, size_t n)
{
    using std::pow;
    assert(a > 1);
    size_t const s = inth_root(n, a);
    auto const mu = mobiusSieve(s);
    auto sieve = range(0, s, [&](size_t k) { return pow(T(k), b) * mu[k]; });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](size_t k) -> T { return sieve[inth_root(k, a)]; }};
}

/// ζ(s) / ζ(2s). f(n) = |μ(n)| = [n is squarefree]. O(n^(3/5)). Motive = -[-1].
template <typename T = u64> Dirichlet<T> squarefree(size_t n, double alpha = 0.6)
{
    size_t const s =
        std::min(Dirichlet<T>::pivot_max, std::max(Dirichlet<T>::defaultPivot(n), (size_t)(alpha * std::pow(n, 0.6))));
    auto const mu = mobiusSieve(isqrt(n));
    auto const M = partialSum(mu, 0);
    auto const precomp = partialSum(std::execution::par, squarefreeSieve(s), T{});
    Dirichlet<T> res(n);
    std::copy_n(std::execution::par, precomp.begin(), res.up().size(), res.up().begin());
    tbb::parallel_for(1UZ, res.down().size(), [&](size_t i) {
        size_t const k = n / i;
        if (k < precomp.size())
            res.down(i) = precomp[k];
        else
        {
            u32 const c_k = inth_root(k, 3);
            res.down(i) = sum(1, c_k, [&](u32 j) { return T(mu[j]) * (k / ((size_t)j * j)) + M[isqrt(k / j)]; }) -
                          T(c_k) * M[c_k];
        }
    });
    return res;
}

/// ζ(s)^2. f(n) = number of divisors of n. O(n^(2/3)). Motive = 2[1].
template <typename T = u64> Dirichlet<T> tau(size_t n, double alpha = 0.08)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivot_max, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    auto precomp = SPF{s}.divisorCountSieve<T>(s);
    partialSumInPlace(precomp);
    Dirichlet<T> S{n};
    std::copy_n(std::execution::par, precomp.begin(), S.up().size(), S.up().begin());
    u32 const max_i = n / precomp.size();
    for (u32 i = max_i + 1; i < S.down().size(); ++i)
        S.down(i) = precomp[S.quotient(i)];
    tbb::parallel_for(1U, max_i + 1, [&](u32 i) {
        libdivide::divider<size_t> const fast_i(i);
        size_t const k = n / i;
        u32 const s = isqrt(k);
        u32 const mid = (S.down().size() - 1) / i;
        S.down(i) = 2 * (sum(1, mid, [&](u32 j) -> T { return S.quotient(i * j); }) +
                         sum(mid + 1, s, [&](u32 j) -> T { return S.quotient(j) / fast_i; })) -
                    T(s) * T(s);
    });
    return S;
}

/// ζ(s)ζ(s - 1). f(n) = sum of divisors of n. O(n^(2/3)). Motive = [p] + [1].
template <typename T = u64> Dirichlet<T> sigma(size_t n, double alpha = 0.08)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivot_max, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    auto precomp = SPF{s}.divisorSumSieve<T>(s);
    partialSumInPlace(precomp);
    Dirichlet<T> S{n};
    std::copy_n(std::execution::par, precomp.begin(), S.up().size(), S.up().begin());
    u32 const max_i = n / precomp.size();
    for (u32 i = max_i + 1; i < S.down().size(); ++i)
        S.down(i) = precomp[S.quotient(i)];
    tbb::parallel_for(1U, max_i + 1, [&](u32 i) {
        libdivide::divider<u64> const fast_i(i);
        u64 const k = n / i;
        u32 const s = isqrt(k);
        u32 const mid = (S.down().size() - 1) / i;
        for (u32 j = 1; j <= mid; ++j)
        {
            u64 const k = S.quotient(i * j);
            S.down(i) += j * k + sumId<T>(k);
        }
        for (u32 j = mid + 1; j <= s; ++j)
        {
            u64 const k = S.quotient(j) / fast_i;
            S.down(i) += j * k + sumId<T>(k);
        }
        S.down(i) -= s * sumId<T>(s);
    });
    return S;
}

/// ζ(s)ζ(s - 2). f(n) = sum of squares of divisors of n. O(n^(2/3)). Motive = [p^2] + [1].
template <typename T = u64> Dirichlet<T> sigma2(size_t n) { return id2<T>(n).multiply(zeta<T>()); }

/// ζ(s)ζ(s - 3). f(n) = sum of cubes of divisors of n. O(n^(2/3)). Motive = [p^3] + [1].
template <typename T = u64> Dirichlet<T> sigma3(size_t n) { return id3<T>(n).multiply(zeta<T>()); }

/// 1 / ζ(s). f(n) = μ(n). F(n) is the Mertens function. O(n^(2/3)). Motive = -[1].
template <typename T = int> Dirichlet<T> mobius(size_t n, double alpha = 0.15)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivot_max, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    auto up = SPF{s}.mobiusSieve<T>(s);
    partialSumInPlace(up);
    return unit<T>(n).divide(zeta<T>(), std::move(up));
}

/// ζ(s - 1) / ζ(s). f(n) = φ(n). O(n^(2/3)). Motive = [p] - [1].
template <typename T = u64> Dirichlet<T> totient(size_t n, double alpha = 0.15)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivot_max, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    auto up = SPF{s}.totientSieve<T>(s);
    partialSumInPlace(up);
    return id<T>(n).divide(zeta<T>(), std::move(up));
}

/// ζ(2s) / ζ(s). f(n) = (-1)^(number of primes dividing n). O(n^(2/3)). Motive = [-1].
template <typename T = int> Dirichlet<T> liouville(size_t n, double alpha = 0.15)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivot_max, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    auto up = SPF{s}.liouvilleSieve<T>(s);
    partialSumInPlace(up);
    return zeta_2s<T>(n).divide(zeta<T>(), std::move(up));
}
} // namespace euler::dirichlet
