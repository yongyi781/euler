#pragma once

#include "decls.hpp"
#include "euler/math.hpp"
#include "it/base.hpp"
#include "libdivide.h"

inline namespace euler
{
namespace dirichlet
{
/// Gives the size of the up vector for a given value of `n`.
constexpr size_t threshold(size_t n)
{
    size_t const s = 0.5 * std::pow(n / log(n), 2.0 / 3);
    // Impose maximum so as to not overwhelm memory.
    size_t const hs = 1'000'000'000;
    return std::min(hs, std::max(isqrt(n), n / (n / s)));
}
} // namespace dirichlet

/// Class for computing Dirichlet series summatory functions, where `up` has size (n / log(n))^(2/3).
template <typename T = int64_t> class Dirichlet
{
  public:
    Dirichlet() = default;

    constexpr Dirichlet(size_t n)
        : _n(n), _up(dirichlet::threshold(n) + 1), _down(n / (dirichlet::threshold(n) + 1) + 1),
          _quotients(isqrt(n) + 1)
    {
        for (size_t i = 1; i < _quotients.size(); ++i)
            _quotients[i] = n / i;
    }

    /// Initialize with the given summatory function.
    template <typename Fun> constexpr Dirichlet(size_t n, Fun F) : Dirichlet(n)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] = F(k);
        for (uint32_t i = _down.size() - 1; i > 0; --i)
            _down[i] = F(_quotients[i]);
    }

    constexpr T &operator[](size_t i) { return i < _up.size() ? _up[i] : _down[_n / i]; }
    constexpr const T &operator[](size_t i) const { return i < _up.size() ? _up[i] : _down[_n / i]; }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] constexpr size_t n() const { return _n; }
    /// The up vector.
    [[nodiscard]] std::vector<T> &up() { return _up; }
    /// The up vector.
    [[nodiscard]] constexpr const std::vector<T> &up() const { return _up; }
    /// The down vector.
    [[nodiscard]] std::vector<T> &down() { return _down; }
    /// The down vector.
    [[nodiscard]] constexpr const std::vector<T> &down() const { return _down; }
    [[nodiscard]] constexpr T &front() { return _up[1]; }
    [[nodiscard]] constexpr const T &front() const { return _up[1]; }
    [[nodiscard]] constexpr T &back() { return _down[1]; }
    [[nodiscard]] constexpr const T &back() const { return _down[1]; }

    /// Enumerates keys of this floors array in ascending order. Breaks if `f` returns `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t ascending(Fun f) const
    {
        for (size_t k = 1; k < _up.size(); ++k)
            if (!it::callbackResult(f, k))
                return it::result_break;
        for (uint32_t i = _down.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, _quotients[i]))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (keys, mutable value) pairs of this floors array in ascending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t ascendingMut(Fun f)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            if (!it::callbackResult(f, k, _up[k]))
                return it::result_break;
        for (uint32_t i = _down.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, _quotients[i], _down[i]))
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
        for (size_t k = _up.size() - 1; k > 0; --k)
            if (!it::callbackResult(f, k))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (key, mutable value) pairs of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t descendingMut(Fun f)
    {
        for (uint32_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, _quotients[i], _down[i]))
                return it::result_break;
        for (size_t k = _up.size() - 1; k > 0; --k)
            if (!it::callbackResult(f, k, _up[k]))
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
        return 2 * (sumMaybeParallel(1, u, [&](uint32_t j) { return (_up[j] - _up[j - 1]) * _down[i * j]; }) +
                    sumMaybeParallel(u + 1, s,
                                     [&](uint32_t j) { return (_up[j] - _up[j - 1]) * _up[_quotients[j] / fasti]; })) -
               _up[s] * _up[s];
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
                                [&](uint32_t j) {
                                    return (_up[j] - _up[j - 1]) * other.down()[i * j] +
                                           (other.up()[j] - other.up()[j - 1]) * _down[i * j];
                                }) +
               sumMaybeParallel(u + 1, s,
                                [&](uint32_t j) {
                                    auto const kdivj = _quotients[j] / fasti;
                                    return (_up[j] - _up[j - 1]) * other.up()[kdivj] +
                                           (other.up()[j] - other.up()[j - 1]) * _up[kdivj];
                                }) -
               _up[s] * other.up()[s];
    }

    constexpr Dirichlet &squareInPlace()
    {
        for (uint32_t i = 1; i < _down.size(); ++i)
            _down[i] = squareValue(_quotients[i]);
        // Sieve for the up values.
        std::vector<T> sieve(_up.size());
        for (size_t i = 1; i < _up.size(); ++i)
        {
            auto const a = _up[i] - _up[i - 1];
            for (size_t j = 1; i * j < _up.size(); ++j)
                sieve[i * j] += a * (_up[j] - _up[j - 1]);
        }
        _up = sieve;
        partialSumInPlace(_up);
        return *this;
    }

    [[nodiscard]] constexpr Dirichlet square() const { return Dirichlet{*this}.squareInPlace(); }

    /// Division in place.
    /// Precondition: `precomputed` must have size at least √n, or be empty.
    template <typename U, std::ranges::sized_range Range = std::ranges::empty_view<T>>
    constexpr Dirichlet &divideInPlace(const Dirichlet<U> &other, Range &&precomputed = {})
    {
        if (precomputed.empty())
            return *this /= other;
        assert(_n == other.n());
        std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
        uint32_t const s = _n / precomputed.size();
        for (uint32_t i = s + 1; i < _down.size(); ++i)
            _down[i] = precomputed[_quotients[i]];
        for (uint32_t i = s; i != 0; --i)
            quotientStep(other, i);
        return *this;
    }

    /// Division.
    /// Precondition: `precomputed` must have size at least √n, or be empty.
    template <typename U, std::ranges::sized_range Range = std::ranges::empty_view<T>>
    constexpr Dirichlet divide(const Dirichlet<U> &other, Range &&precomputed = {}) const
    {
        return Dirichlet{*this}.divideInPlace(other, std::forward<Range>(precomputed));
    }

    /// Addition.
    template <typename U> constexpr Dirichlet &operator+=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] += other.up()[k];
        for (uint32_t i = 1; i < _down.size(); ++i)
            _down[i] += other.down()[i];
        return *this;
    }
    template <typename U> constexpr Dirichlet operator+(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} += other;
    }

    /// Subtraction.
    template <typename U> constexpr Dirichlet &operator-=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] -= other.up()[k];
        for (uint32_t i = 1; i < _down.size(); ++i)
            _down[i] -= other.down()[i];
        return *this;
    }
    template <typename U> constexpr Dirichlet operator-(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} -= other;
    }

    /// Multiplication.
    template <typename U> constexpr Dirichlet &operator*=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (uint32_t i = 1; i < _down.size(); ++i)
            _down[i] = productValue(other, _n / i);
        // Sieve for the up values.
        std::vector<T> sieve(_up.size());
        for (size_t i = 1; i < _up.size(); ++i)
        {
            auto const a = _up[i] - _up[i - 1];
            for (size_t j = 1; i * j < _up.size(); ++j)
                sieve[i * j] += a * (other.up()[j] - other.up()[j - 1]);
        }
        _up = sieve;
        partialSumInPlace(_up);
        return *this;
    }
    template <typename U> constexpr Dirichlet operator*(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} *= other;
    }

    /// Division.
    template <typename U> constexpr Dirichlet &operator/=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        auto const c = other.up()[1] - other.up()[0];
        // Sieve for the up values.
        adjacentDifferenceInPlace(_up);
        for (size_t i = 1; i <= (_up.size() - 1) / 2; ++i)
        {
            auto const a = _up[i];
            for (size_t j = 2; i * j < _up.size(); ++j)
                _up[i * j] -= a * (other.up()[j] - other.up()[j - 1]);
            _up[i] /= c;
        }
        partialSumInPlace(_up);
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            quotientStep(other, i);
        return *this;
    }
    template <typename U> constexpr Dirichlet operator/(const Dirichlet<U> &other) const
    {
        return Dirichlet{*this} /= other;
    }

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
        _down[i] -= (_up[1] - _up[0]) * other[k];
        _down[i] -=
            sumMaybeParallel(2, u,
                             [&](uint32_t j) {
                                 return (_up[j] - _up[j - 1]) * other.down()[i * j] +
                                        (other.up()[j] - other.up()[j - 1]) * _down[i * j];
                             }) +
            sumMaybeParallel(u + 1, s, [&](uint32_t j) {
                size_t const kdivj = _quotients[j] / fasti;
                return (_up[j] - _up[j - 1]) * other.up()[kdivj] + (other.up()[j] - other.up()[j - 1]) * _up[kdivj];
            });
        _down[i] += _up[s] * other.up()[s];
        _down[i] /= other.up()[1] - other.up()[0];
    }
};

namespace dirichlet
{
/// 1, the multiplicative identity. f(n) = [n = 1].
template <typename T = int64_t> constexpr Dirichlet<T> unit(size_t n)
{
    return {n, [](auto &&) { return T(1); }};
}

/// ζ(s). f(n) = 1.
template <typename T = int64_t> constexpr Dirichlet<T> zeta(size_t n) { return {n, std::identity{}}; }

/// ζ(s - 1). f(n) = n.
template <typename T = int64_t> constexpr Dirichlet<T> id(size_t n)
{
    return {n, [](auto &&k) { return T(k) * (k + 1) / 2; }};
}

/// ζ(s - 2). f(n) = n^2.
template <typename T = int64_t> constexpr Dirichlet<T> id2(size_t n)
{
    return {n, [](auto &&k) { return T(k) * (k + 1) * (2 * k + 1) / 6; }};
}

/// ζ(2s). f(n) = [n is square].
template <typename T = int64_t> constexpr Dirichlet<T> zeta_2s(size_t n)
{
    return {n, [](auto &&k) { return isqrt(k); }};
}

/// ζ(rs). f(n) = [n is a perfect rth power].
template <typename T = int64_t> constexpr Dirichlet<T> zeta_rs(size_t n, int r)
{
    return {n, [&](auto &&k) { return (T)std::pow(k + 0.5, 1.0 / r); }};
}

/// χ_4(s). f(n) = (-4|n). 1 if 1 mod 4, -1 if 3 mod 4, 0 otherwise. Also known as the Dirichlet beta function.
template <typename T = int64_t> constexpr Dirichlet<T> chi4(size_t n)
{
    return {n, [&](auto &&k) { return T(k % 4 == 1 || k % 4 == 2); }};
}

/// 1 / ζ(s). f(n) = μ(n). F(n) is the Mertens function. O(n^(2/3)).
template <typename T = int64_t> constexpr Dirichlet<T> mu(size_t n)
{
    size_t const s = std::max(threshold(n), (size_t)(0.35 * std::pow(n, 2.0 / 3)));
    auto sieve = mobiusSieve(s);
    auto ssieve = partialSum(sieve, T{});
    auto res = unit<T>(n);
    res.divideInPlace(zeta<T>(n), ssieve);
    return res;
}

/// ζ(s) / ζ(2s). f(n) = |μ(n)| = [n is squarefree]. O(n^(2/3)).
template <typename T = int64_t> constexpr Dirichlet<T> squarefree(size_t n)
{
    size_t const s = std::max(threshold(n), (size_t)(0.6 * std::pow(n, 2.0 / 3)));
    auto sieve = squarefreeSieve(s);
    auto ssieve = partialSum(sieve, T{});
    auto res = zeta<T>(n);
    res.divideInPlace(zeta_2s<T>(n), ssieve);
    return res;
}

/// ζ(s - 1) / ζ(s). f(n) = φ(n). O(n^(2/3)).
template <typename T = int64_t> constexpr Dirichlet<T> totient(size_t n)
{
    size_t const s = std::max(threshold(n), (size_t)(0.6 * std::pow(n, 2.0 / 3)));
    auto sieve = totientSieve(s);
    auto ssieve = partialSum(sieve, T{});
    auto res = id<T>(n);
    res.divideInPlace(zeta<T>(n), ssieve);
    return res;
}

/// ζ(s)^2. f(n) = number of divisors of n.
template <typename T = int64_t> constexpr Dirichlet<T> tau(size_t n)
{
    // auto sieve = divisorCountSieve((size_t)(0.85 * std::pow(n, 2.0 / 3)));
    // auto ssieve = partialSum(sieve, T{});
    auto res = zeta<T>(n);
    res.squareInPlace();
    return res;
}

/// ζ(s)ζ(s-1). f(n) = sum of divisors of n.
template <typename T = int64_t> constexpr Dirichlet<T> sigma(size_t n)
{
    // auto sieve = divisorCountSieve((size_t)(0.85 * std::pow(n, 2.0 / 3)));
    // auto ssieve = partialSum(sieve, T{});
    return zeta<T>(n) * id<T>(n);
}

/// ζ(2s) / ζ(s). f(n) = (-1)^(number of primes dividing n).
template <typename T = int64_t> constexpr Dirichlet<T> liouville(size_t n)
{
    // auto sieve = divisorCountSieve((size_t)(0.85 * std::pow(n, 2.0 / 3)));
    // auto ssieve = partialSum(sieve, T{});
    return zeta_2s<T>(n) / zeta<T>(n);
}
} // namespace dirichlet
} // namespace euler
