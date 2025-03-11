#pragma once

#include "decls.hpp"
#include "euler/math.hpp"
#include "it/base.hpp"
#include "libdivide.h"

inline namespace euler
{
/// Class for computing Dirichlet series summatory functions, where `up` has size O((n / log(n))^(2/3)).
template <typename T = int64_t> class Dirichlet
{
  public:
    Dirichlet() = default;

    /// Configurable upper bound on the size of the up vector.
    static inline size_t pivotMax = 1'000'000'000;
    /// Configurable coefficient on the size of the up vector.
    static inline double pivotCoefficient = 0.5;

    /// Gives the default size of the up vector for a given value of `n`.
    static constexpr size_t defaultPivot(size_t n)
    {
        size_t const s = pivotCoefficient * std::pow(n / log(n), 2.0 / 3);
        // Impose maximum so as to not overwhelm memory.
        size_t const hs = pivotMax;
        size_t const res = std::max(isqrt(n), std::min(hs, s));
        return n / (n / res);
    }

    constexpr explicit Dirichlet(size_t n)
        : _n(n), _up(defaultPivot(n) + 1), _down(n / (defaultPivot(n) + 1) + 1), _quotients(isqrt(n) + 1)
    {
        for (size_t i = 1; i < _quotients.size(); ++i)
            _quotients[i] = n / i;
    }

    /// Initialize with the given summatory function.
    template <typename Fun> constexpr Dirichlet(size_t n, Fun F) : Dirichlet(n)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] = F(k);
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            _down[i] = F(_quotients[i]);
    }

    constexpr T &operator[](size_t i) { return i < _up.size() ? _up[i] : _down[_n / i]; }
    constexpr const T &operator[](size_t i) const { return i < _up.size() ? _up[i] : _down[_n / i]; }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] constexpr size_t n() const { return _n; }
    /// The index of the transition point between up and down vectors.
    [[nodiscard]] constexpr size_t pivot() const { return _up.size() - 1; }
    /// Gets the value of `n / i` avoiding a division CPU instruction. `i` must be `≤ √n`.
    [[nodiscard]] constexpr size_t quotient(uint32_t i) const { return _quotients[i]; }

    /// The up vector, mutable.
    [[nodiscard]] constexpr std::vector<T> &up() { return _up; }
    /// The up vector.
    [[nodiscard]] constexpr const std::vector<T> &up() const { return _up; }
    /// Element access into the up array, mutable.
    [[nodiscard]] constexpr T &up(size_t x) { return _up[x]; }
    /// Element access into the up array.
    [[nodiscard]] constexpr const T &up(size_t x) const { return _up[x]; }

    /// The down vector, mutable.
    [[nodiscard]] constexpr std::vector<T> &down() { return _down; }
    /// The down vector.
    [[nodiscard]] constexpr const std::vector<T> &down() const { return _down; }
    /// Element access into the down array, mutable.
    [[nodiscard]] constexpr T &down(size_t x) { return _down[x]; }
    /// Element access into the down array.
    [[nodiscard]] constexpr const T &down(size_t x) const { return _down[x]; }

    [[nodiscard]] constexpr T &front() { return up(1); }
    [[nodiscard]] constexpr const T &front() const { return up(1); }
    [[nodiscard]] constexpr T &back() { return down(1); }
    [[nodiscard]] constexpr const T &back() const { return down(1); }

    /// Returns a vector of function values from 1 to the size of the up vector. This is the vector of adjacent
    /// differences of `up`.
    [[nodiscard]] constexpr std::vector<T> upValues() const { return adjacentDifference(_up); }

    /// Takes a partial sum in place. Useful function to call when building this array from individual values.
    constexpr Dirichlet &accumulate()
    {
        for (size_t i = 1; i < _up.size(); ++i)
            up(i) += up(i - 1);
        _down.back() += _up.back();
        for (uint32_t i = _down.size() - 2; i != 0; --i)
            down(i) += down(i + 1);
        return *this;
    }

    /// Applies a function to each entry of this array.
    template <std::invocable<T> Fun> constexpr Dirichlet &mapInPlace(Fun f)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) = f(up(k));
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            down(i) = f(down(i));
        return *this;
    }

    /// Enumerates keys of this floors array in ascending order. Breaks if `f` returns `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t ascending(Fun f) const
    {
        for (size_t k = 1; k < _up.size(); ++k)
            if (!it::callbackResult(f, k))
                return it::result_break;
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            if (!it::callbackResult(f, _quotients[i]))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (keys, mutable value) pairs of this floors array in ascending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t ascendingMut(Fun f)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            if (!it::callbackResult(f, k, up(k)))
                return it::result_break;
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            if (!it::callbackResult(f, _quotients[i], down(i)))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates keys of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t descending(Fun f) const
    {
        for (uint32_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, _quotients[i]))
                return it::result_break;
        for (size_t k = _up.size() - 1; k != 0; --k)
            if (!it::callbackResult(f, k))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (key, mutable value) pairs of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t descendingMut(Fun f)
    {
        for (uint32_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, _quotients[i], down(i)))
                return it::result_break;
        for (size_t k = _up.size() - 1; k != 0; --k)
            if (!it::callbackResult(f, k, up(k)))
                return it::result_break;
        return it::result_continue;
    }

    /// Compute a single summatory value of `this * this`.
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    [[nodiscard]] constexpr T squareValue(size_t k) const
    {
        uint32_t const s = isqrt(k);
        uint32_t const i = _n / k;
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return 2 * (sumMaybeParallel(1, u, [&](uint32_t j) -> T { return (up(j) - up(j - 1)) * down(i * j); }) +
                    sumMaybeParallel(
                        u + 1, s, [&](uint32_t j) -> T { return (up(j) - up(j - 1)) * up(_quotients[j] / fasti); })) -
               up(s) * up(s);
    }

    /// Compute a single summatory value of (this * other).
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <typename U> [[nodiscard]] constexpr T productValue(const Dirichlet<U> &other, size_t k) const
    {
        if constexpr (std::is_same_v<T, U>)
            if (this == &other)
                return squareValue(k);
        uint32_t const s = isqrt(k);
        uint32_t const i = _n / k;
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return sumMaybeParallel(1, u,
                                [&](uint32_t j) -> T {
                                    return (up(j) - up(j - 1)) * other.down(i * j) +
                                           (other.up(j) - other.up(j - 1)) * down(i * j);
                                }) +
               sumMaybeParallel(u + 1, s,
                                [&](uint32_t j) -> T {
                                    auto const kdivj = _quotients[j] / fasti;
                                    return (up(j) - up(j - 1)) * other.up(kdivj) +
                                           (other.up(j) - other.up(j - 1)) * up(kdivj);
                                }) -
               up(s) * other.up(s);
    }

    constexpr Dirichlet &squareInPlace()
    {
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) = squareValue(_quotients[i]);
        // Sieve for the up values.
        std::vector<T> sieve(_up.size());
        uint32_t const s = isqrt(_up.size() - 1);
        T a{};
        for (size_t i = 1; i <= s; ++i)
        {
            a = up(i) - up(i - 1);
            sieve[i * i] += a * a;
            for (size_t j = i + 1; i * j < _up.size(); ++j)
                sieve[i * j] += 2 * a * (up(j) - up(j - 1));
        }
        _up = sieve;
        partialSumInPlace(_up);
        return *this;
    }

    [[nodiscard]] constexpr Dirichlet square() const { return Dirichlet{*this}.squareInPlace(); }

    /// Division in place.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <typename U, std::ranges::sized_range Range = std::ranges::empty_view<T>>
    constexpr Dirichlet &divideInPlace(const Dirichlet<U> &other, Range &&precomputed = {})
    {
        if (precomputed.empty())
            return *this /= other;
        assert(_n == other.n());
        std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
        uint32_t const s = _n / precomputed.size();
        for (uint32_t i = s + 1; i < _down.size(); ++i)
            down(i) = precomputed[_quotients[i]];
        for (uint32_t i = s; i != 0; --i)
            quotientStep(other, i);
        return *this;
    }

    /// Division.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <typename U, std::ranges::sized_range Range = std::ranges::empty_view<T>>
    [[nodiscard]] constexpr Dirichlet divide(const Dirichlet<U> &other, Range &&precomputed = {}) const
    {
        return Dirichlet{*this}.divideInPlace(other, std::forward<Range>(precomputed));
    }

    /// Dirichlet inverse.
    [[nodiscard]] constexpr Dirichlet inverse() const
    {
        Dirichlet res{_n, [](auto &&) { return T(1); }};
        return res /= *this;
    }

    /// Addition.
    template <typename U> constexpr Dirichlet &operator+=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) += other.up(k);
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) += other.down(i);
        return *this;
    }
    template <typename U> [[nodiscard]] constexpr Dirichlet operator+(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} += other;
    }

    /// Subtraction.
    template <typename U> constexpr Dirichlet &operator-=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) -= other.up(k);
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) -= other.down(i);
        return *this;
    }
    template <typename U> [[nodiscard]] constexpr Dirichlet operator-(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} -= other;
    }

    /// Multiplication.
    template <typename U> constexpr Dirichlet &operator*=(const Dirichlet<U> &other)
    {
        if constexpr (std::is_same_v<T, U>)
            if (this == &other)
                return squareInPlace();
        assert(_n == other.n());
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) = productValue(other, _n / i);
        // Sieve for the up values.
        std::vector<T> sieve(_up.size());
        T a{};
        for (size_t i = 1; i < _up.size(); ++i)
        {
            a = up(i) - up(i - 1);
            for (size_t j = 1; i * j < _up.size(); ++j)
                sieve[i * j] += a * (other.up(j) - other.up(j - 1));
        }
        _up = sieve;
        partialSumInPlace(_up);
        return *this;
    }
    template <typename U> [[nodiscard]] constexpr Dirichlet operator*(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} *= other;
    }

    /// Division.
    template <typename U> constexpr Dirichlet &operator/=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        auto const c = other.up(1) - other.up(0);
        // Sieve for the up values.
        adjacentDifferenceInPlace(_up);
        for (size_t i = 1; i <= (_up.size() - 1) / 2; ++i)
        {
            if (c != 1)
                up(i) /= c;
            for (size_t j = 2; i * j < _up.size(); ++j)
                up(i * j) -= up(i) * (other.up(j) - other.up(j - 1));
        }
        if (c != 1)
            for (size_t k = (_up.size() + 1) / 2; k < _up.size(); ++k)
                up(k) /= c;
        partialSumInPlace(_up);
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            quotientStep(other, i);
        return *this;
    }
    template <typename U> [[nodiscard]] constexpr Dirichlet operator/(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} /= other;
    }

    /// Dirichlet inverse.
    [[nodiscard]] constexpr Dirichlet operator~() const { return inverse(); }

    template <typename U> constexpr bool operator==(const Dirichlet<U> &other) const
    {
        return _n == other.n() && _up == other.up() && _down == other.down();
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Dirichlet &S)
    {
        return o << "{\n  n: " << S._n << "\n  up: " << S._up << "\n  down: " << S._down << "\n}";
    }

  private:
    size_t _n = 0;
    std::vector<T> _up;
    std::vector<T> _down;
    std::vector<size_t> _quotients;

    /// One step of the divison algorithm. Internal use only.
    template <typename U> constexpr void quotientStep(const Dirichlet<U> &other, uint32_t i)
    {
        size_t const k = _n / i;
        uint32_t const s = isqrt(k);
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        down(i) -= (up(1) - up(0)) * other[k];
        down(i) -= sumMaybeParallel(2, u,
                                    [&](uint32_t j) -> T {
                                        return (up(j) - up(j - 1)) * other.down(i * j) +
                                               (other.up(j) - other.up(j - 1)) * down(i * j);
                                    }) +
                   sumMaybeParallel(u + 1, s, [&](uint32_t j) -> T {
                       size_t const kdivj = _quotients[j] / fasti;
                       return (up(j) - up(j - 1)) * other.up(kdivj) + (other.up(j) - other.up(j - 1)) * up(kdivj);
                   });
        down(i) += up(s) * other.up(s);
        down(i) /= other.up(1) - other.up(0);
    }
};

namespace dirichlet
{
/// 1, the multiplicative identity. f(n) = [n = 1]. Motive = 0.
template <typename T = int64_t> constexpr Dirichlet<T> unit(size_t n)
{
    return {n, [](auto &&) { return T(1); }};
}

/// ζ(s). f(n) = 1. Motive = [1].
template <typename T = int64_t> constexpr Dirichlet<T> zeta(size_t n) { return {n, std::identity{}}; }

/// ζ(s - 1). f(n) = n. Motive = [p].
template <typename T = int64_t> constexpr Dirichlet<T> id(size_t n)
{
    return {n, [](auto &&k) { return k % 2 == 0 ? T(k / 2) * (k + 1) : T((k + 1) / 2) * k; }};
}

/// ζ(s - 2). f(n) = n^2. If `T = ZMod<M>`, currently requires that neither 2 or 3 divide `M`. Motive = [p^2].
template <typename T = int64_t> constexpr Dirichlet<T> id2(size_t n)
{
    return {n, [](auto &&k) { return T(k) * (k + 1) * (2 * k + 1) / 6; }};
}

/// ζ(s - 3). f(n) = n^3. Motive = [p^3].
template <typename T = int64_t> constexpr Dirichlet<T> id3(size_t n)
{
    return {n, [](auto &&k) {
                T const s = k % 2 == 0 ? T(k / 2) * (k + 1) : T((k + 1) / 2) * k;
                return s * s;
            }};
}

/// ζ(2s). f(n) = [n is square]. Motive = [1] + [-1].
template <typename T = int64_t> constexpr Dirichlet<T> zeta_2s(size_t n)
{
    return {n, [](auto &&k) { return isqrt(k); }};
}

/// ζ(as). f(n) = [n is a perfect ath power]. Motive = sum([r] for r in ath roots of unity)
template <typename T = int64_t> constexpr Dirichlet<T> zeta_multiple(size_t n, int a)
{
    return {n, [&](auto &&k) { return (T)std::pow(k + 0.5, 1.0 / a); }};
}

/// χ_(-3)(s). f(n) = (-3|n). 1 if 1 mod 3, -1 if 2 mod 3, 0 otherwise. L-function of the hexagonal lattice and the
/// Eisenstein integers.
template <typename T = int64_t> constexpr Dirichlet<T> chi3(size_t n)
{
    return {n, [&](auto &&k) { return T(k % 3 == 1); }};
}

/// χ_(-4)(s). f(n) = (-4|n). 1 if 1 mod 4, -1 if 3 mod 4, 0 otherwise. Also known as the Dirichlet beta function.
/// Motive = χ_(-4).
template <typename T = int64_t> constexpr Dirichlet<T> chi4(size_t n)
{
    return {n, [&](auto &&k) { return T(k % 4 == 1 || k % 4 == 2); }};
}

/// χ_5(s). f(n) = (5|n).
template <typename T = int64_t> constexpr Dirichlet<T> chi5(size_t n)
{
    return {n, [&](auto &&k) { return std::array{0, 1, 0, -1, 0}[k % 5]; }};
}

/// 1 / ζ(2s). f(n) = [n is square] * μ(√n). O(n^(1/2)). Motive = -[1] - [-1].
template <typename T = int64_t> constexpr Dirichlet<T> inv_zeta_2s(size_t n)
{
    auto const mu = mobiusSieve(isqrt(n));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](auto &&k) { return mertens[isqrt(k)]; }};
}

/// 1 / ζ(as). f(n) = [n is a perfect rth power] * μ(n^(1/a)). O(n^(1/a)). Motive = -∑([r] for r in ath roots of unity).
template <typename T = int64_t> constexpr Dirichlet<T> inv_zeta_multiple(size_t n, int a)
{
    auto const mu = mobiusSieve((size_t)std::pow(n + 0.5, 1.0 / a));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](auto &&k) { return mertens[(size_t)std::pow(k + 0.5, 1.0 / a)]; }};
}

/// ζ(as - b). f(n) = [n = k^a] * k^b. Requires `a > 1`. Motive = ∑([r*p^(b/a)] for r in ath roots of unity).
template <typename T = int64_t> constexpr Dirichlet<T> zeta_linear(size_t n, int a, int b)
{
    assert(a > 1);
    size_t const s = std::pow(n + 0.5, 1.0 / a);
    auto sieve = range(0, s, [&](size_t k) { return std::pow(T(k), b); });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](auto &&k) { return sieve[std::pow(k + 0.5, 1.0 / a)]; }};
}

/// 1 / ζ(as - b). f(n) = [n = k^a] * μ(k) * k^b. Requires `a > 1`. Motive = -∑([r*p^(b/a)] for r in ath roots of
/// unity).
template <typename T = int64_t> constexpr Dirichlet<T> inv_zeta_linear(size_t n, int a, int b)
{
    assert(a > 1);
    size_t const s = std::pow(n + 0.5, 1.0 / a);
    auto const mu = mobiusSieve((size_t)std::pow(n + 0.5, 1.0 / a));
    auto sieve = range(0, s, [&](size_t k) { return std::pow(T(k), b) * mu[k]; });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](auto &&k) { return sieve[std::pow(k + 0.5, 1.0 / a)]; }};
}

/// ζ(s) / ζ(2s). f(n) = |μ(n)| = [n is squarefree]. O(n^(3/5)). Motive = -[-1].
template <typename T = int64_t> constexpr Dirichlet<T> squarefree(size_t n)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(1.25 * std::pow(n, 0.6)));
    auto const mu = mobiusSieve(isqrt(n));
    auto const mertens = partialSum(mu, T{});
    auto const precomputed = partialSum(squarefreeSieve(s), T{});
    return {n, [&](size_t k) -> T {
                if (k < precomputed.size())
                    return precomputed[k];
                uint32_t const s = cbrt(k);
                return sum(1, s,
                           [&](uint32_t j) { return mu[j] * (k / ((size_t)j * j)) + mertens[std::sqrt(k / j)]; }) -
                       mertens[s] * s;
            }};
}

/// 1 / ζ(s). f(n) = μ(n). F(n) is the Mertens function. O(n^(2/3)). Motive = -[1].
template <typename T = int64_t> constexpr Dirichlet<T> mu(size_t n)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(0.35 * std::pow(n, 2.0 / 3)));
    return unit<T>(n).divideInPlace(zeta<T>(n), partialSum(mobiusSieve(s), T{}));
}

/// ζ(s - 1) / ζ(s). f(n) = φ(n). O(n^(2/3)). Motive = [p] - [1].
template <typename T = int64_t> constexpr Dirichlet<T> totient(size_t n)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(0.6 * std::pow(n, 2.0 / 3)));
    return id<T>(n).divideInPlace(zeta<T>(n), partialSum(totientSieve(s), T{}));
}

/// ζ(s)^2. f(n) = number of divisors of n. O(n^(2/3)). Motive = 2[1].
template <typename T = int64_t> constexpr Dirichlet<T> tau(size_t n) { return zeta<T>(n).squareInPlace(); }

/// ζ(s)ζ(s - 1). f(n) = sum of divisors of n. O(n^(2/3)). Motive = [p] + [1].
template <typename T = int64_t> constexpr Dirichlet<T> sigma(size_t n) { return zeta<T>(n) * id<T>(n); }

/// ζ(s)ζ(s - 2). f(n) = sum of squares of divisors of n. O(n^(2/3)). Motive = [p^2] + [1].
template <typename T = int64_t> constexpr Dirichlet<T> sigma2(size_t n) { return zeta<T>(n) * id2<T>(n); }

/// ζ(s)ζ(s - 3). f(n) = sum of cubes of divisors of n. O(n^(2/3)). Motive = [p^3] + [1].
template <typename T = int64_t> constexpr Dirichlet<T> sigma3(size_t n) { return zeta<T>(n) * id3<T>(n); }

/// ζ(2s) / ζ(s). f(n) = (-1)^(number of primes dividing n). O(n^(2/3)). Motive = [-1].
template <typename T = int64_t> constexpr Dirichlet<T> liouville(size_t n) { return zeta_2s<T>(n) / zeta<T>(n); }
} // namespace dirichlet
} // namespace euler
