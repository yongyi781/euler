#pragma once

#include "SPF.hpp"
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
    static inline size_t pivotMax = 1UZ << 31;
    /// Configurable exponent on the size of the up vector.
    static inline double pivotExponent = 2.0 / 3;
    /// Configurable coefficient on the size of the up vector. Ideal is 0.2 for multiplication and 0.5 for division.
    static inline double pivotCoefficient = 0.2;

    /// Gives the default size of the up vector for a given value of `n`.
    static constexpr size_t defaultPivot(size_t n)
    {
        size_t const s = pivotCoefficient * std::pow(n / std::max(1.0, log(n)), pivotExponent);
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
    [[nodiscard]] constexpr size_t n() const noexcept { return _n; }
    /// The index of the transition point between up and down vectors.
    [[nodiscard]] constexpr size_t pivot() const noexcept { return _up.size() - 1; }
    /// Gets the value of `n / i` avoiding a division CPU instruction. `i` must be `≤ √n`.
    [[nodiscard]] constexpr size_t quotient(uint32_t i) const noexcept { return _quotients[i]; }

    /// The up vector, mutable.
    [[nodiscard]] constexpr std::vector<T> &up() noexcept { return _up; }
    /// The up vector.
    [[nodiscard]] constexpr const std::vector<T> &up() const noexcept { return _up; }
    /// Element access into the up array, mutable.
    [[nodiscard]] constexpr T &up(size_t x) noexcept { return _up[x]; }
    /// Element access into the up array.
    [[nodiscard]] constexpr const T &up(size_t x) const noexcept { return _up[x]; }

    /// The down vector, mutable.
    [[nodiscard]] constexpr std::vector<T> &down() noexcept { return _down; }
    /// The down vector.
    [[nodiscard]] constexpr const std::vector<T> &down() const noexcept { return _down; }
    /// Element access into the down array, mutable.
    [[nodiscard]] constexpr T &down(size_t x) noexcept { return _down[x]; }
    /// Element access into the down array.
    [[nodiscard]] constexpr const T &down(size_t x) const noexcept { return _down[x]; }

    [[nodiscard]] constexpr T &front() noexcept { return up(1); }
    [[nodiscard]] constexpr const T &front() const noexcept { return up(1); }
    [[nodiscard]] constexpr T &back() noexcept { return _down.size() > 1 ? down(1) : _up.back(); }
    [[nodiscard]] constexpr const T &back() const noexcept { return _down.size() > 1 ? down(1) : _up.back(); }

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
            if (!it::callbackResult(f, quotient(i)))
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
            if (!it::callbackResult(f, quotient(i), down(i)))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates keys of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t descending(Fun f) const
    {
        for (uint32_t i = 1; i < _down.size(); ++i)
            if (!it::callbackResult(f, quotient(i)))
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
            if (!it::callbackResult(f, quotient(i), down(i)))
                return it::result_break;
        for (size_t k = _up.size() - 1; k != 0; --k)
            if (!it::callbackResult(f, k, up(k)))
                return it::result_break;
        return it::result_continue;
    }

    /// Compute a single summatory value of `this * this`.
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <size_t ParThreshold = 8192> [[nodiscard]] constexpr T squareValue(size_t k) const
    {
        uint32_t const s = isqrt(k);
        uint32_t const i = _n / k;
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return 2 * (sumMaybeParallel<ParThreshold>(1, u,
                                                   [&](uint32_t j) -> T { return (up(j) - up(j - 1)) * down(i * j); }) +
                    sumMaybeParallel<ParThreshold>(
                        u + 1, s, [&](uint32_t j) -> T { return (up(j) - up(j - 1)) * up(quotient(j) / fasti); })) -
               up(s) * up(s);
    }

    /// Compute a single summatory value of (this * other).
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <size_t ParThreshold = 8192, typename U>
    [[nodiscard]] constexpr T productValue(const Dirichlet<U> &other, size_t k) const
    {
        if constexpr (std::is_same_v<T, U>)
            if (this == &other)
                return squareValue(k);
        uint32_t const s = isqrt(k);
        uint32_t const i = _n / k;
        uint32_t const u = (_down.size() - 1) / i;
        if (u < 1 || s < u + 1)
            std::cout << "Help";
        libdivide::divider<size_t> const fasti(i);
        return sumMaybeParallel<ParThreshold>(1, u,
                                              [&](uint32_t j) -> T {
                                                  return (up(j) - up(j - 1)) * other.down(i * j) +
                                                         (other.up(j) - other.up(j - 1)) * down(i * j);
                                              }) +
               sumMaybeParallel<ParThreshold>(u + 1, s,
                                              [&](uint32_t j) -> T {
                                                  auto const kdivj = quotient(j) / fasti;
                                                  return (up(j) - up(j - 1)) * other.up(kdivj) +
                                                         (other.up(j) - other.up(j - 1)) * up(kdivj);
                                              }) -
               up(s) * other.up(s);
    }

    /// Compute a single summatory value of (this * ζ). About 33% faster than generic.
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <size_t ParThreshold = 8192> [[nodiscard]] constexpr T productValueZeta(size_t k) const
    {
        uint32_t const s = isqrt(k);
        uint32_t const i = _n / k;
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return sumMaybeParallel<ParThreshold>(
                   1, u, [&](uint32_t j) -> T { return (up(j) - up(j - 1)) * quotient(i * j) + down(i * j); }) +
               sumMaybeParallel<ParThreshold>(u + 1, s,
                                              [&](uint32_t j) -> T {
                                                  auto const kdivj = quotient(j) / fasti;
                                                  return (up(j) - up(j - 1)) * kdivj + up(kdivj);
                                              }) -
               up(s) * s;
    }

    /// Squares this Dirichlet series in place.
    template <std::ranges::sized_range Range = std::ranges::empty_view<T>>
    constexpr Dirichlet &squareInPlace(Range &&precomputed = {})
    {
        uint32_t const u = precomputed.empty() ? _down.size() - 1 : _n / precomputed.size();
        _down = mapv(std::execution::par, range(0_u32, (uint32_t)_down.size() - 1),
                     [&](uint32_t i) { return i == 0 ? T(0) : squareValue<0>(quotient(i)); });

        if (precomputed.empty())
        {
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
        }
        else
        {
            for (uint32_t i = u + 1; i < _down.size(); ++i)
                down(i) = precomputed[quotient(i)];
            std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
        }
        return *this;
    }

    /// Multiplies this Dirichlet series with ζ, in place. About 33% faster than generic.
    template <std::ranges::sized_range Range = std::ranges::empty_view<T>>
    constexpr Dirichlet &multiplyInPlaceZeta(Range &&precomputed = {})
    {
        uint32_t const u = precomputed.empty() ? _down.size() - 1 : _n / precomputed.size();
        _down = mapv(std::execution::par, range(0_u32, (uint32_t)_down.size() - 1),
                     [&](uint32_t i) { return i == 0 ? T(0) : productValueZeta<0>(quotient(i)); });
        if (precomputed.empty())
        {
            // Sieve for the up values.
            std::vector<T> sieve(_up.size());
            T a{};
            for (size_t i = 1; i < _up.size(); ++i)
            {
                a = up(i) - up(i - 1);
                for (size_t j = 1; i * j < _up.size(); ++j)
                    sieve[i * j] += a;
            }
            _up = sieve;
            partialSumInPlace(_up);
        }
        else
        {
            for (uint32_t i = u + 1; i < _down.size(); ++i)
                down(i) = precomputed[quotient(i)];
            std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
        }
        return *this;
    }

    [[nodiscard]] constexpr Dirichlet square(this Dirichlet self)
    {
        self.squareInPlace();
        return self;
    }

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
            down(i) = precomputed[quotient(i)];
        for (uint32_t i = s; i != 0; --i)
            quotientStep(other, i);
        return *this;
    }

    /// Division in place, in the special case that the divisor is ζ(s). 33% speedup compared to generic.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <std::ranges::sized_range Range = std::ranges::empty_view<T>>
    constexpr Dirichlet &divideInPlaceZeta(Range &&precomputed = {})
    {
        uint32_t u = _down.size() - 1;
        if (precomputed.empty())
        {
            // Sieve for the up values.
            adjacentDifferenceInPlace(_up);
            for (size_t i = 1; i <= (_up.size() - 1) / 2; ++i)
                for (size_t j = 2; i * j < _up.size(); ++j)
                    up(i * j) -= up(i);
            partialSumInPlace(_up);
        }
        else
        {
            std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
            u = _n / precomputed.size();
            for (uint32_t i = u + 1; i < _down.size(); ++i)
                down(i) = precomputed[quotient(i)];
        }
        for (uint32_t i = u; i != 0; --i)
            quotientStepZeta(i);
        return *this;
    }

    /// Division.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <typename U, std::ranges::sized_range Range = std::ranges::empty_view<T>>
    [[nodiscard]] constexpr Dirichlet divide(this Dirichlet self, const Dirichlet<U> &other, Range &&precomputed = {})
    {
        self.divideInPlace(other, std::forward<Range>(precomputed));
        return self;
    }

    /// Dirichlet inverse.
    [[nodiscard]] constexpr Dirichlet inverse() const
    {
        Dirichlet res{_n, [](auto &&) { return T(1); }};
        res /= *this;
        return res;
    }

    /// Returns this raised to a power.
    [[nodiscard]] constexpr Dirichlet pow(int exponent) const
    {
        Dirichlet x{_n, [](auto &&) { return T(1); }};
        if (exponent == 0)
            return x;
        Dirichlet y = *this;
        if (exponent < 0)
        {
            y = ~y;
            exponent = -exponent;
        }
        while (true)
        {
            if (exponent & 1)
            {
                x *= y;
                if (exponent == 1)
                    break;
            }
            exponent >>= 1;
            y.squareInPlace();
        }
        return x;
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

    template <typename U> [[nodiscard]] constexpr friend Dirichlet operator+(Dirichlet left, const Dirichlet<U> &right)
    {
        left += right;
        return left;
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

    template <typename U> [[nodiscard]] constexpr friend Dirichlet operator-(Dirichlet left, const Dirichlet<U> &right)
    {
        left -= right;
        return left;
    }

    /// Multiplication.
    template <typename U> constexpr Dirichlet &operator*=(const Dirichlet<U> &other)
    {
        if constexpr (std::is_same_v<T, U>)
            if (this == &other)
                return squareInPlace();
        assert(_n == other.n());
        _down = mapv(std::execution::par, range(0_u32, (uint32_t)_down.size() - 1),
                     [&](uint32_t i) { return i == 0 ? T(0) : productValue<0>(other, quotient(i)); });
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

    template <typename U> [[nodiscard]] constexpr friend Dirichlet operator*(Dirichlet left, const Dirichlet<U> &right)
    {
        left *= right;
        return left;
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

    template <typename U> [[nodiscard]] constexpr friend Dirichlet operator/(Dirichlet left, const Dirichlet<U> &right)
    {
        left /= right;
        return left;
    }

    /// Multiplication by a scalar.
    constexpr Dirichlet &operator*=(T value)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) *= value;
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) *= value;
        return *this;
    }

    [[nodiscard]] constexpr friend Dirichlet operator*(Dirichlet left, T value)
    {
        left *= value;
        return left;
    }

    /// Division by a scalar.
    constexpr Dirichlet &operator/=(T value)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) /= value;
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) /= value;
        return *this;
    }

    [[nodiscard]] constexpr friend Dirichlet operator/(Dirichlet left, T value)
    {
        left /= value;
        return left;
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
        size_t const k = quotient(i);
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
                       size_t const kdivj = quotient(j) / fasti;
                       return (up(j) - up(j - 1)) * other.up(kdivj) + (other.up(j) - other.up(j - 1)) * up(kdivj);
                   });
        down(i) += up(s) * other.up(s);
        down(i) /= other.up(1) - other.up(0);
    }

    /// One step of the divison algorithm by ζ(s). Internal use only.
    constexpr void quotientStepZeta(uint32_t i)
    {
        size_t const k = quotient(i);
        uint32_t const s = isqrt(k);
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        down(i) -= (up(1) - up(0)) * k;
        down(i) -= sumMaybeParallel(
                       2, u, [&](uint32_t j) -> T { return (up(j) - up(j - 1)) * quotient(i * j) + down(i * j); }) +
                   sumMaybeParallel(u + 1, s, [&](uint32_t j) -> T {
                       size_t const kdivj = quotient(j) / fasti;
                       return (up(j) - up(j - 1)) * kdivj + up(kdivj);
                   });
        down(i) += up(s) * s;
    }
};

namespace dirichlet
{
/// Computes `∑ (ab ≤ n), g(a) * h(b)` in O(√n).
/// Also known as `∑ (k ≤ n), (g * h)(k)`.
/// Also known as `∑ (k ≤ n), g(k) * H(n/k) = ∑ (k ≤ n), h(k) * G(n/k)`.
/// Requirements:
/// * `G` and `H` should be the summatory functions of `g` and `h`.
/// * Need to be able to evaluate `g(k)`, `h(k)` for `k ≤ √n` and `G(m)` and `H(m)` for `m ≥ √n`.
template <integral2 T = size_t, typename Fun1, typename SummatoryFun1, typename Fun2, typename SummatoryFun2>
auto productValue(T n, Fun1 f, SummatoryFun1 F, Fun2 g, SummatoryFun2 G)
{
    using H = half_integer_t<T>;
    using Tp = std::common_type_t<std::invoke_result_t<SummatoryFun1, T>, std::invoke_result_t<SummatoryFun1, T>>;
    H const s = isqrt(n);
    return sumMaybeParallel(H(1), s,
                            [&](H k) {
                                T const ndivk = n / k;
                                return Tp(F(ndivk)) * g(k) + Tp(G(ndivk)) * f(k);
                            }) -
           Tp(F(s)) * G(s);
}

/// Computes `∑ (ab ≤ n), g(a) * h(b)` in O(√n), given only their summatory functions. Convenience function.
template <integral2 Tn = size_t, typename SummatoryFun1, typename SummatoryFun2>
auto productValue(Tn n, SummatoryFun1 F, SummatoryFun2 G)
{
    return productValue(
        std::move(n), [&](auto &&k) { return F(k) - F(k - 1); }, F, [&](auto &&k) { return G(k) - G(k - 1); }, G);
}

/// 1, the multiplicative identity. f(n) = [n = 1]. Motive = 0.
template <typename T = int64_t> constexpr Dirichlet<T> unit(size_t n)
{
    return {n, [](size_t) { return T(1); }};
}

/// ζ(s). f(n) = 1. Motive = [1].
template <typename T = int64_t> constexpr Dirichlet<T> zeta(size_t n) { return {n, std::identity{}}; }

/// ζ(s - 1). f(n) = n. Motive = [p].
template <typename T = int64_t> constexpr Dirichlet<T> id(size_t n)
{
    return {n, [](size_t k) { return sumId<T>(k); }};
}

/// ζ(s - 2). f(n) = n^2. If `T = ZMod<M>`, currently requires that neither 2 or 3 divide `M`. Motive = [p^2].
template <typename T = int64_t> constexpr Dirichlet<T> id2(size_t n)
{
    // Fancy stuff to avoid overflow or ZMod divisions.
    return {n, [](size_t k) { return sumSquares<T>(k); }};
}

/// ζ(s - 3). f(n) = n^3. Motive = [p^3].
template <typename T = int64_t> constexpr Dirichlet<T> id3(size_t n)
{
    return {n, [](size_t k) {
                T const s = k % 2 == 0 ? T(k / 2) * (k + 1) : T((k + 1) / 2) * k;
                return s * s;
            }};
}

/// ζ(2s). f(n) = [n is square]. Motive = [1] + [-1].
template <typename T = int64_t> constexpr Dirichlet<T> zeta_2s(size_t n)
{
    return {n, [](size_t k) { return isqrt(k); }};
}

/// ζ(as). f(n) = [n is a perfect ath power]. Motive = sum([r] for r in ath roots of unity)
template <typename T = int64_t> constexpr Dirichlet<T> zeta_multiple(int a, size_t n)
{
    return {n, [&](size_t k) { return (T)std::pow(k + 0.5, 1.0 / a); }};
}

/// χ_(-3)(s). f(n) = (-3|n). 1 if 1 mod 3, -1 if 2 mod 3, 0 otherwise. L-function of the hexagonal lattice and the
/// Eisenstein integers.
template <typename T = int64_t> constexpr Dirichlet<T> chi3(size_t n)
{
    return {n, [&](size_t k) { return T(k % 3 == 1); }};
}

/// χ_(-4)(s). f(n) = (-4|n). 1 if 1 mod 4, -1 if 3 mod 4, 0 otherwise. Also known as the Dirichlet beta function.
/// Motive = χ_(-4).
template <typename T = int64_t> constexpr Dirichlet<T> chi4(size_t n)
{
    return {n, [&](size_t k) { return T(k % 4 == 1 || k % 4 == 2); }};
}

/// χ_5(s). f(n) = (5|n).
template <typename T = int64_t> constexpr Dirichlet<T> chi5(size_t n)
{
    return {n, [&](size_t k) { return std::array{0, 1, 0, -1, 0}[k % 5]; }};
}

/// 1 / ζ(2s). f(n) = [n is square] * μ(√n). O(n^(1/2)). Motive = -[1] - [-1].
template <typename T = int64_t> constexpr Dirichlet<T> inv_zeta_2s(size_t n)
{
    auto const mu = mobiusSieve(isqrt(n));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](size_t k) { return mertens[isqrt(k)]; }};
}

/// 1 / ζ(as). f(n) = [n is a perfect rth power] * μ(n^(1/a)). O(n^(1/a)). Motive = -∑([r] for r in ath roots of unity).
template <typename T = int64_t> constexpr Dirichlet<T> inv_zeta_multiple(int a, size_t n)
{
    auto const mu = mobiusSieve((size_t)std::pow(n + 0.5, 1.0 / a));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](size_t k) { return mertens[(size_t)std::pow(k + 0.5, 1.0 / a)]; }};
}

/// ζ(as - b). f(n) = [n = k^a] * k^b. Requires `a > 1`. Motive = ∑([r*p^(b/a)] for r in ath roots of unity).
template <typename T = int64_t> constexpr Dirichlet<T> zeta_linear(int a, int b, size_t n)
{
    assert(a > 1);
    size_t const s = std::pow(n + 0.5, 1.0 / a);
    auto sieve = range(0, s, [&](size_t k) { return std::pow(T(k), b); });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](size_t k) { return sieve[std::pow(k + 0.5, 1.0 / a)]; }};
}

/// 1 / ζ(as - b). f(n) = [n = k^a] * μ(k) * k^b. Requires `a > 1`. Motive = -∑([r*p^(b/a)] for r in ath roots of
/// unity).
template <typename T = int64_t> constexpr Dirichlet<T> inv_zeta_linear(int a, int b, size_t n)
{
    assert(a > 1);
    size_t const s = std::pow(n + 0.5, 1.0 / a);
    auto const mu = mobiusSieve((size_t)std::pow(n + 0.5, 1.0 / a));
    auto sieve = range(0, s, [&](size_t k) { return std::pow(T(k), b) * mu[k]; });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](size_t k) { return sieve[std::pow(k + 0.5, 1.0 / a)]; }};
}

/// ζ(s) / ζ(2s). f(n) = |μ(n)| = [n is squarefree]. O(n^(3/5)). Motive = -[-1].
template <typename T = int64_t> constexpr Dirichlet<T> squarefree(size_t n, double alpha = 1.25)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(alpha * std::pow(n, 0.6)));
    auto const mu = mobiusSieve(isqrt(n));
    auto const mertens = partialSum(mu, T{});
    auto const precomputed = partialSum(squarefreeSieve<uint8_t>(s), T{});
    return {n, [&](size_t k) -> T {
                if (k < precomputed.size())
                    return precomputed[k];
                uint32_t const s = cbrt(k);
                return sum(1, s,
                           [&](uint32_t j) { return mu[j] * (k / ((size_t)j * j)) + mertens[std::sqrt(k / j)]; }) -
                       mertens[s] * s;
            }};
}

/// ζ(s)^2. f(n) = number of divisors of n. O(n^(2/3)). Motive = 2[1].
template <typename T = int64_t> constexpr Dirichlet<T> tau(size_t n, double alpha = 0.05)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(alpha * std::pow(n, 2.0 / 3)));
    auto const precomputed = partialSum(divisorCountSieve(s), T{});
    Dirichlet S{n};
    std::copy(precomputed.begin(), precomputed.begin() + S.up().size(), S.up().begin());
    uint32_t const u = n / precomputed.size();
    for (uint32_t i = u + 1; i < S.down().size(); ++i)
        S.down(i) = precomputed[S.quotient(i)];
    std::for_each(std::execution::par, counting_iterator(1_u32), counting_iterator(u + 1), [&](uint32_t i) {
        size_t const k = n / i;
        uint32_t const s = isqrt(k);
        uint32_t const u = (S.down().size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        S.down(i) = 2 * (sum(1, u, [&](uint32_t j) -> T { return S.quotient(i * j); }) +
                         sum(u + 1, s, [&](uint32_t j) -> T { return S.quotient(j) / fasti; })) -
                    T(s) * T(s);
    });
    return S;
}

/// ζ(s)ζ(s - 1). f(n) = sum of divisors of n. O(n^(2/3)). Motive = [p] + [1].
template <typename T = int64_t> constexpr Dirichlet<T> sigma(size_t n, double alpha = 0.38)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(alpha * std::pow(n, 2.0 / 3)));
    return id<T>(n).multiplyInPlaceZeta(partialSum(divisorSumSieve(s), T{}));
}

/// ζ(s)ζ(s - 2). f(n) = sum of squares of divisors of n. O(n^(2/3)). Motive = [p^2] + [1].
template <typename T = int64_t> constexpr Dirichlet<T> sigma2(size_t n) { return id2<T>(n).multiplyInPlaceZeta(); }

/// ζ(s)ζ(s - 3). f(n) = sum of cubes of divisors of n. O(n^(2/3)). Motive = [p^3] + [1].
template <typename T = int64_t> constexpr Dirichlet<T> sigma3(size_t n) { return id3<T>(n).multiplyInPlaceZeta(); }

/// 1 / ζ(s). f(n) = μ(n). F(n) is the Mertens function. O(n^(2/3)). Motive = -[1].
template <typename T = int64_t> constexpr Dirichlet<T> mobius(size_t n, double alpha = 0.38)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(alpha * std::pow(n, 2.0 / 3)));
    return unit<T>(n).divideInPlaceZeta(partialSum(mobiusSieve(s), T{}));
}

/// ζ(s - 1) / ζ(s). f(n) = φ(n). O(n^(2/3)). Motive = [p] - [1].
template <typename T = int64_t> constexpr Dirichlet<T> totient(size_t n, double alpha = 0.5)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(alpha * std::pow(n, 2.0 / 3)));
    return id<T>(n).divideInPlaceZeta(partialSum(totientSieve(s), T{}));
}

/// ζ(2s) / ζ(s). f(n) = (-1)^(number of primes dividing n). O(n^(2/3)). Motive = [-1].
template <typename T = int64_t> constexpr Dirichlet<T> liouville(size_t n, double alpha = 0.5)
{
    size_t const s = std::max(Dirichlet<>::defaultPivot(n), (size_t)(alpha * std::pow(n, 2.0 / 3)));
    return zeta_2s<T>(n).divideInPlaceZeta(partialSum(liouvilleSieve(s), T{}));
}
} // namespace dirichlet
} // namespace euler
