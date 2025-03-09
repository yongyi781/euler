#pragma once

#include "decls.hpp"
#include "euler/math.hpp"
#include "it/base.hpp"

inline namespace euler
{
/// Class for computing Dirichlet series summatory functions.
template <typename T = int64_t> class Dirichlet
{
  public:
    Dirichlet() = default;
    constexpr Dirichlet(size_t n) : _n(n), _up(isqrt(n) + 1), _down(n / (isqrt(n) + 1) + 1) {}
    /// Initialize with the given summatory function.
    template <typename Fun> constexpr Dirichlet(size_t n, Fun F) : Dirichlet(n)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] = F(k);
        for (size_t i = _down.size() - 1; i > 0; --i)
            _down[i] = F(_n / i);
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
        for (size_t i = 1; i < _up.size(); ++i)
            if (!it::callbackResult(f, i))
                return it::result_break;
        for (size_t i = _down.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, _n / i))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (keys, mutable value) pairs of this floors array in ascending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t ascendingMut(Fun f)
    {
        for (size_t i = 1; i < _up.size(); ++i)
            if (!it::callbackResult(f, i, _up[i]))
                return it::result_break;
        for (size_t i = _down.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, _n / i, _down[i]))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates keys of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t descending(Fun f) const
    {
        for (size_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, _n / i))
                return it::result_break;
        for (size_t i = _up.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, i))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (key, mutable value) pairs of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t descendingMut(Fun f)
    {
        for (size_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, _n / i, _down[i]))
                return it::result_break;
        for (size_t i = _up.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, i, _up[i]))
                return it::result_break;
        return it::result_continue;
    }

    /// Compute a single summatory value of (this * this).
    [[nodiscard]] constexpr T squareValue(size_t n) const
    {
        size_t const s = isqrt(n);
        auto res = 2 * sum(1, s, [&](auto &&k) { return (_up[k] - _up[k - 1]) * (*this)[n / k]; });
        res -= _up[s] * _up[s];
        return res;
    }

    /// Compute a single summatory value of (this * other).
    template <typename U> [[nodiscard]] constexpr T productValue(const Dirichlet<U> &other, size_t n) const
    {
        if constexpr (std::is_same_v<T, U>)
            if (this == &other)
                return squareValue(n);
        size_t const s = isqrt(n);
        auto res = sumMaybeParallel<8000>(1, s, [&](auto &&k) {
            auto const ndivk = n / k;
            return (_up[k] - _up[k - 1]) * other[ndivk] + (*this)[ndivk] * (other.up()[k] - other.up()[k - 1]);
        });
        res -= _up[s] * other.up()[s];
        return res;
    }

    constexpr Dirichlet &squareInPlace()
    {
        for (size_t i = 1; i < _down.size(); ++i)
        {
            size_t const k = _n / i;
            size_t const s = isqrt(k);
            size_t const u = (_down.size() - 1) / i;
            _down[i] = 2 * (sum(1, u, [&](auto &&j) { return (_up[j] - _up[j - 1]) * _down[i * j]; }) +
                            sum(u + 1, s, [&](auto &&j) { return (_up[j] - _up[j - 1]) * _up[k / j]; })) -
                       _up[s] * _up[s];
        }
        for (size_t k = _up.size() - 1; k > 0; --k)
        {
            size_t const s = isqrt(k);
            _up[k] = 2 * sum(1, s, [&](auto &&j) { return (_up[j] - _up[j - 1]) * _up[k / j]; }) - _up[s] * _up[s];
        }
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
        auto const c = other[1] - other[0];
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] = precomputed[k];
        size_t const u = _n / precomputed.size();
        for (size_t i = u + 1; i < _down.size(); ++i)
            _down[i] = precomputed[_n / i];
        for (size_t i = u; i != 0; --i)
        {
            size_t const k = _n / i;
            size_t const s = isqrt(k);
            size_t const v = (_down.size() - 1) / i;
            _down[i] -=
                sumMaybeParallel<10000>(2, v,
                                        [&](auto &&j) {
                                            return (_up[j] - _up[j - 1]) * other.down()[i * j] +
                                                   (other.up()[j] - other.up()[j - 1]) * _down[i * j];
                                        }) +
                sumMaybeParallel<10000>(v + 1, s, [&](auto &&j) {
                    auto const kdivj = k / j;
                    return (_up[j] - _up[j - 1]) * other.up()[kdivj] + (other.up()[j] - other.up()[j - 1]) * _up[kdivj];
                });

            _down[i] -= (_up[1] - _up[0]) * other[k];
            _down[i] += _up[s] * other.up()[s];
            _down[i] /= c;
        }
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
        for (size_t i = 1; i < _down.size(); ++i)
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
        for (size_t i = 1; i < _down.size(); ++i)
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
        for (size_t i = 1; i < _down.size(); ++i)
        {
            size_t const k = _n / i;
            size_t const s = isqrt(k);
            size_t const u = (_down.size() - 1) / i;
            _down[i] = (sumMaybeParallel<10000>(1, u,
                                                [&](auto &&j) {
                                                    return (_up[j] - _up[j - 1]) * other.down()[i * j] +
                                                           (other.up()[j] - other.up()[j - 1]) * _down[i * j];
                                                }) +
                        sumMaybeParallel<10000>(u + 1, s,
                                                [&](auto &&j) {
                                                    auto const kdivj = k / j;
                                                    return (_up[j] - _up[j - 1]) * other.up()[kdivj] +
                                                           (other.up()[j] - other.up()[j - 1]) * _up[kdivj];
                                                })) -
                       _up[s] * other.up()[s];
        }
        for (size_t k = _up.size() - 1; k > 0; --k)
        {
            size_t const s = isqrt(k);
            _up[k] = sumMaybeParallel<10000>(1, s,
                                             [&](auto &&j) {
                                                 auto const kdivj = k / j;
                                                 return (_up[j] - _up[j - 1]) * other.up()[kdivj] +
                                                        (other.up()[j] - other.up()[j - 1]) * _up[kdivj];
                                             }) -
                     _up[s] * other.up()[s];
        }
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
        auto const c = other[1] - other[0];
        _up[1] /= c;
        for (size_t k = 2; k < _up.size(); ++k)
        {
            auto const s = isqrt(k);
            _up[k] -= (_up[1] - _up[0]) * other[k];
            _up[k] -= sumMaybeParallel<10000>(2, s, [&](auto &&j) {
                auto const kdivj = k / j;
                return (_up[j] - _up[j - 1]) * other.up()[kdivj] + _up[kdivj] * (other.up()[j] - other.up()[j - 1]);
            });
            _up[k] += _up[s] * other.up()[s];
            _up[k] /= c;
        }
        for (size_t i = _down.size() - 1; i != 0; --i)
        {
            auto const k = _n / i;
            auto const s = isqrt(k);
            _down[i] -= (_up[1] - _up[0]) * other[k];
            _down[i] -= sumMaybeParallel<10000>(2, s, [&](auto &&j) {
                auto const kdivj = k / j;
                return (_up[j] - _up[j - 1]) * other[kdivj] + (*this)[kdivj] * (other.up()[j] - other.up()[j - 1]);
            });
            _down[i] += _up[s] * other.up()[s];
            _down[i] /= c;
        }
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

/// 1 / ζ(s). f(n) = μ(n). O(n^(2/3)).
template <typename T = int64_t> constexpr Dirichlet<T> mu(size_t n)
{
    auto sieve = mobiusSieve((size_t)(0.85 * std::pow(n, 2.0 / 3)));
    auto ssieve = partialSum(sieve, T{});
    auto res = unit<T>(n);
    res.divideInPlace(zeta<T>(n), ssieve);
    return res;
}

/// ζ(s) / ζ(2s). f(n) = [n is squarefree]. O(n^(2/3)).
template <typename T = int64_t> constexpr Dirichlet<T> squarefree(size_t n)
{
    auto sieve = squarefreeSieve(std::pow(n, 2.0 / 3));
    auto ssieve = partialSum(sieve, T{});
    auto res = zeta<T>(n);
    res.divideInPlace(zeta_2s<T>(n), ssieve);
    return res;
}

/// ζ(s - 1) / ζ(s). f(n) = φ(n). O(n^(2/3)).
template <typename T = int64_t> constexpr Dirichlet<T> totient(size_t n)
{
    auto sieve = totientSieve((size_t)(0.85 * std::pow(n, 2.0 / 3)));
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

/// ζ(2s) / ζ(s). f(n) = sum of divisors of n.
template <typename T = int64_t> constexpr Dirichlet<T> liouville(size_t n)
{
    // auto sieve = divisorCountSieve((size_t)(0.85 * std::pow(n, 2.0 / 3)));
    // auto ssieve = partialSum(sieve, T{});
    return zeta_2s<T>(n) / zeta<T>(n);
}
} // namespace dirichlet
} // namespace euler
