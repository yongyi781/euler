#pragma once

#include "SPF.hpp"
#include "decls.hpp"
#include "euler/math.hpp"
#include "it/tree.hpp"
#include "libdivide.h"
#include "literals.hpp"

namespace euler
{
template <typename T = int64_t> class Dirichlet;
template <typename Fun, typename SFun> class SpecialDirichlet;

template <typename T> struct is_dirichlet : std::false_type
{
};

template <typename T> struct is_dirichlet<Dirichlet<T>> : std::true_type
{
};

template <typename Fun, typename SFun> struct is_dirichlet<SpecialDirichlet<Fun, SFun>> : std::true_type
{
};

template <typename T>
concept dirichlet_type = is_dirichlet<std::decay_t<T>>::value;

/// A class representing a special Dirichlet series that does not require any storage.
template <typename Fun, typename SFun> class SpecialDirichlet
{
    Fun m_f;
    SFun m_F;

  public:
    constexpr SpecialDirichlet(Fun f, SFun F) : m_f(std::move(f)), m_F(std::move(F)) {}

    [[nodiscard]] constexpr auto value(size_t k) const noexcept { return m_f(k); }
    [[nodiscard]] constexpr auto sum(size_t k) const noexcept { return m_F(k); }
    [[nodiscard]] constexpr auto up(size_t k) const noexcept { return m_F(k); }

    [[nodiscard]] constexpr auto operator[](size_t k) const noexcept { return m_F(k); }
    [[nodiscard]] constexpr auto operator()(size_t k) const noexcept { return m_F(k); }

    template <typename Fun2, typename SFun2, integral2 T>
    [[nodiscard]] auto productValue(const SpecialDirichlet<Fun2, SFun2> &other, T n) const
    {
        using H = half_integer_t<T>;
        using Tp = std::common_type_t<std::invoke_result_t<SFun, T>, std::invoke_result_t<SFun2, T>>;
        H const s = isqrt(n);
        return sumMaybeParallel(H(1), s,
                                [&](H k) {
                                    T const ndivk = n / k;
                                    return Tp(sum(ndivk)) * Tp(other.value(k)) + Tp(other.sum(ndivk)) * Tp(value(k));
                                }) -
               Tp(sum(s)) * other.sum(s);
    }
};

/// Class for computing Dirichlet series summatory functions, where `up` has size O((n / log(n))^(2/3)).
template <typename T> class Dirichlet
{
    size_t _n = 0;
    std::vector<T> _up;
    std::vector<T> _down;
    std::vector<size_t> _quotients;

    template <dirichlet_type Dir> T sumTerm1(const Dir &other, size_t i, uint32_t j) const
    {
        if constexpr (requires { other.down(i * j); })
            return value(j) * other.down(i * j) + down(i * j) * other.value(j);
        return value(j) * other.sum(quotient(i * j)) + down(i * j) * other.value(j);
    }

    template <dirichlet_type Dir> T sumTerm2(const Dir &other, const libdivide::divider<size_t> &i, uint32_t j) const
    {
        size_t const kdivj = quotient(j) / i;
        return value(j) * other.up(kdivj) + up(kdivj) * other.value(j);
    }

    /// One step of the divison algorithm by ζ(s). Internal use only.
    template <dirichlet_type Dir> void divideStep(const Dir &other, uint32_t i)
    {
        size_t const k = quotient(i);
        uint32_t const s = isqrt(k);
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        down(i) -= (up(1) - up(0)) * T(other.sum(k));
        down(i) -= sumMaybeParallel(2, u, [&](uint32_t j) -> T { return sumTerm1(other, i, j); });
        down(i) -= sumMaybeParallel(u + 1, s, [&](uint32_t j) -> T { return sumTerm2(other, fasti, j); });
        down(i) += up(s) * T(other.up(s));
        if (other.value(1) != 1)
            down(i) /= other.value(1);
    }

  public:
    Dirichlet() = default;

    /// Default max memory usage = 16 GB.
    static constexpr size_t DefaultMaxMemoryUsage = 1UZ << 34;
    /// Configurable upper bound on the size of the up vector.
    static inline size_t pivotMax = DefaultMaxMemoryUsage / sizeof(T);
    /// Configurable exponent on the size of the up vector.
    static inline double pivotExponent = 2.0 / 3;
    /// Configurable coefficient on the size of the up vector. Ideal is 0.2 for multiplication and 0.5 for division.
    static inline double pivotCoefficient = 0.2;

    /// Gives the default size of the up vector for a given value of `n`.
    static size_t defaultPivot(size_t n)
    {
        size_t const s = pivotCoefficient * std::pow(n / std::max(1.0, log(n)), pivotExponent);
        // Impose maximum so as to not overwhelm memory.
        size_t const hs = pivotMax;
        size_t const res = std::max(isqrt(n), std::min(hs, s));
        return n / (n / res);
    }

    explicit Dirichlet(size_t n)
        : _n(n), _up(defaultPivot(n) + 1), _down(n / (defaultPivot(n) + 1) + 1), _quotients(isqrt(n) + 1)
    {
        for (size_t i = 1; i < _quotients.size(); ++i)
            _quotients[i] = n / i;
    }

    /// Initialize with the given summatory function.
    template <typename SummatoryFun> Dirichlet(size_t n, SummatoryFun F) : Dirichlet(n)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            _up[k] = T(F(k));
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            _down[i] = T(F(_quotients[i]));
    }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] size_t n() const noexcept { return _n; }
    /// The index of the transition point between up and down vectors.
    [[nodiscard]] size_t pivot() const noexcept { return _up.size() - 1; }
    /// Gets the value of `n / i` avoiding a division CPU instruction. `i` must be `≤ √n`.
    [[nodiscard]] size_t quotient(uint32_t i) const noexcept { return _quotients[i]; }

    /// The up vector, mutable.
    [[nodiscard]] std::vector<T> &up() noexcept { return _up; }
    /// The up vector.
    [[nodiscard]] const std::vector<T> &up() const noexcept { return _up; }
    /// Element access into the up array, mutable.
    [[nodiscard]] T &up(size_t x) noexcept { return _up[x]; }
    /// Element access into the up array.
    [[nodiscard]] const T &up(size_t x) const noexcept { return _up[x]; }

    /// The down vector, mutable.
    [[nodiscard]] std::vector<T> &down() noexcept { return _down; }
    /// The down vector.
    [[nodiscard]] const std::vector<T> &down() const noexcept { return _down; }
    /// Element access into the down array, mutable.
    [[nodiscard]] T &down(size_t x) noexcept { return _down[x]; }
    /// Element access into the down array.
    [[nodiscard]] const T &down(size_t x) const noexcept { return _down[x]; }

    [[nodiscard]] T &front() noexcept { return up(1); }
    [[nodiscard]] const T &front() const noexcept { return up(1); }
    [[nodiscard]] T &back() noexcept { return _down.size() > 1 ? down(1) : _up.back(); }
    [[nodiscard]] const T &back() const noexcept { return _down.size() > 1 ? down(1) : _up.back(); }

    /// Returns a vector of function values from 1 to the pivot. This is the vector of adjacent differences of `up`.
    [[nodiscard]] std::vector<T> values() const { return adjacentDifference(_up); }
    /// Returns a single value. Requires `k ≤ pivot()`.
    [[nodiscard]] T value(size_t k) const noexcept { return up(k) - up(k - 1); }

    [[nodiscard]] T &operator[](size_t k) noexcept { return k < _up.size() ? _up[k] : _down[_n / k]; }
    [[nodiscard]] const T &operator[](size_t k) const noexcept { return k < _up.size() ? _up[k] : _down[_n / k]; }
    [[nodiscard]] T &operator()(size_t k) noexcept { return (*this)[k]; }
    [[nodiscard]] const T &operator()(size_t k) const noexcept { return (*this)[k]; }

    /// Returns the prefix sum of this series at the given index.
    [[nodiscard]] const T &sum(size_t k) const noexcept { return (*this)[k]; }

    /// Takes a partial sum in place. Useful function to call when building this array from individual values.
    Dirichlet &accumulate()
    {
        for (size_t i = 1; i < _up.size(); ++i)
            up(i) += up(i - 1);
        _down.back() += _up.back();
        for (uint32_t i = _down.size() - 2; i != 0; --i)
            down(i) += down(i + 1);
        return *this;
    }

    /// Applies a function to each entry of this array.
    template <std::invocable<T> Fun> Dirichlet &mapInPlace(Fun f)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) = f(up(k));
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            down(i) = f(down(i));
        return *this;
    }

    /// Enumerates keys of this floors array in ascending order. Breaks if `f` returns `it::result_break`.
    template <std::invocable<size_t> Fun> it::result_t ascending(Fun f) const
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
    template <std::invocable<size_t, T &> Fun> it::result_t ascendingMut(Fun f)
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
    template <std::invocable<size_t> Fun> it::result_t descending(Fun f) const
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
    template <std::invocable<size_t, T &> Fun> it::result_t descendingMut(Fun f)
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
    template <size_t ParThreshold = 8192> [[nodiscard]] T squareValue(size_t k) const
    {
        uint32_t const s = isqrt(k);
        size_t const i = _n / k;
        uint32_t const u = (_down.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return 2 * (sumMaybeParallel<ParThreshold>(1, u, [&](uint32_t j) -> T { return value(j) * down(i * j); }) +
                    sumMaybeParallel<ParThreshold>(
                        u + 1, s, [&](uint32_t j) -> T { return value(j) * up(quotient(j) / fasti); })) -
               up(s) * up(s);
    }

    /// Compute a single summatory value of (this * other), where `other` is a `SpecialDirichlet`.
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <size_t ParThreshold = 8192, dirichlet_type Dir>
    [[nodiscard]] T productValue(const Dir &other, size_t k) const
    {
        if constexpr (std::is_same_v<Dir, Dirichlet>)
            if (this == &other)
                return squareValue(k);
        uint32_t const s = isqrt(k);
        size_t const i = _n / k;
        uint32_t const u = (down().size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return sumMaybeParallel<ParThreshold>(1, u, [&](uint32_t j) -> T { return sumTerm1(other, i, j); }) +
               sumMaybeParallel<ParThreshold>(u + 1, s, [&](uint32_t j) -> T { return sumTerm2(other, fasti, j); }) -
               up(s) * other.up(s);
    }

    /// Computes `∑ (p^e * b ≤ n), f(e) * self(b)`.
    template <typename Fun> [[nodiscard]] T localProductValue(size_t p, Fun f, size_t k) const
    {
        assert(p >= 2);
        T res = 0;
        for (int e = 0; k != 0; ++e, k /= p)
        {
            auto const c = f(e);
            if (c != 0)
                res += c * (*this)[k];
        }
        return res;
    }

    /// Squares this Dirichlet series in place.
    template <std::ranges::range Range = std::ranges::empty_view<T>> Dirichlet &squareInPlace(Range &&precomputed = {})
    {
        uint32_t const u = precomputed.empty() ? _down.size() - 1 : _n / precomputed.size();
        _down = mapv(std::execution::par, range(0_u32, (uint32_t)_down.size() - 1),
                     [&](uint32_t i) { return i == 0 || i > u ? T(0) : squareValue<0>(quotient(i)); });

        if (precomputed.empty())
        {
            // Sieve for the up values.
            auto const old = values();
            std::ranges::fill(_up, T(0));
            uint32_t const s = isqrt(_up.size() - 1);
            for (size_t i = 1; i <= s; ++i)
            {
                up(i * i) += old[i] * old[i];
                for (size_t j = i + 1; i * j < _up.size(); ++j)
                    up(i * j) += 2 * old[i] * old[j];
            }
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

    template <std::ranges::range Range = std::ranges::empty_view<T>>
    [[nodiscard]] Dirichlet square(this Dirichlet self, Range &&precomputed = {})
    {
        self.squareInPlace(std::forward<Range>(precomputed));
        return self;
    }

    /// Multiplication.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    Dirichlet &multiplyInPlace(const Dir &other, Range &&precomputed = {})
    {
        if (precomputed.empty())
            return *this *= other;
        uint32_t const s = _n / precomputed.size();
        _down = mapv(std::execution::par, range(0_u32, (uint32_t)_down.size() - 1), [&](uint32_t i) -> T {
            if (i == 0)
                return T(0);
            if (i > s)
                return precomputed[quotient(i)];
            return productValue<0>(other, quotient(i));
        });
        std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
        return *this;
    }

    /// Division.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    [[nodiscard]] Dirichlet multiply(this Dirichlet self, const Dir &other, Range &&precomputed = {})
    {
        self.multiplyInPlace(other, std::forward<Range>(precomputed));
        return self;
    }

    /// Division in place.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    Dirichlet &divideInPlace(const Dir &other, Range &&precomputed = {})
    {
        if (precomputed.empty())
            return *this /= other;
        std::copy(precomputed.begin(), precomputed.begin() + _up.size(), _up.begin());
        uint32_t const s = _n / precomputed.size();
        for (uint32_t i = s + 1; i < _down.size(); ++i)
            down(i) = precomputed[quotient(i)];
        for (uint32_t i = s; i != 0; --i)
            divideStep(other, i);
        return *this;
    }

    /// Division.
    /// Precondition: `precomputed` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    [[nodiscard]] Dirichlet divide(this Dirichlet self, const Dir &other, Range &&precomputed = {})
    {
        self.divideInPlace(other, std::forward<Range>(precomputed));
        return self;
    }

    /// Dirichlet inverse.
    [[nodiscard]] Dirichlet inverse() const
    {
        Dirichlet res{_n, [](auto &&) { return T(1); }};
        res /= *this;
        return res;
    }

    /// Returns this raised to a power.
    template <integral2 E> [[nodiscard]] Dirichlet pow(this Dirichlet self, E exponent)
    {
        Dirichlet x{self.n(), [](auto &&) { return T(1); }};
        if (exponent == 0)
            return x;
        if (exponent < 0)
        {
            self = ~self;
            exponent = -exponent;
        }
        while (true)
        {
            if (exponent & 1)
            {
                x *= self;
                if (exponent == 1)
                    break;
            }
            exponent >>= 1;
            self.squareInPlace();
        }
        return x;
    }

    /// Addition.
    template <typename U> Dirichlet &operator+=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) += other.up(k);
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) += other.down(i);
        return *this;
    }

    template <typename U> [[nodiscard]] friend Dirichlet operator+(Dirichlet left, const Dirichlet<U> &right)
    {
        left += right;
        return left;
    }

    /// Subtraction.
    template <typename U> Dirichlet &operator-=(const Dirichlet<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) -= other.up(k);
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) -= other.down(i);
        return *this;
    }

    template <typename U> [[nodiscard]] friend Dirichlet operator-(Dirichlet left, const Dirichlet<U> &right)
    {
        left -= right;
        return left;
    }

    /// Multiplication.
    template <dirichlet_type Dir> Dirichlet &operator*=(const Dir &other)
    {
        if constexpr (std::is_same_v<Dir, Dirichlet>)
            if (this == &other)
                return squareInPlace();
        _down = mapv(std::execution::par, range(0_u32, (uint32_t)_down.size() - 1),
                     [&](uint32_t i) { return i == 0 ? T(0) : productValue<0>(other, quotient(i)); });
        // Sieve for the up values.
        T const c = other.value(1);
        adjacentDifferenceInPlace(_up);
        T a{};
        for (size_t i = _up.size() - 1; i != 0; --i)
        {
            a = up(i);
            up(i) *= c;
            for (size_t j = 2; i * j < _up.size(); ++j)
                up(i * j) += a * other.value(j);
        }
        partialSumInPlace(_up);
        return *this;
    }

    template <dirichlet_type Dir> [[nodiscard]] friend Dirichlet operator*(Dirichlet left, const Dir &right)
    {
        left *= right;
        return left;
    }

    template <typename Fun, typename SFun>
    [[nodiscard]] friend Dirichlet operator*(const SpecialDirichlet<Fun, SFun> &left, Dirichlet right)
    {
        right *= left;
        return right;
    }

    /// Division.
    template <dirichlet_type Dir> Dirichlet &operator/=(const Dir &other)
    {
        auto const c = other.value(1);
        // Sieve for the up values.
        adjacentDifferenceInPlace(_up);
        for (size_t i = 1; i <= (_up.size() - 1) / 2; ++i)
        {
            if (c != 1)
                up(i) /= c;
            for (size_t j = 2; i * j < _up.size(); ++j)
                up(i * j) -= up(i) * other.value(j);
        }
        if (c != 1)
            for (size_t k = (_up.size() + 1) / 2; k < _up.size(); ++k)
                up(k) /= c;
        partialSumInPlace(_up);
        for (uint32_t i = _down.size() - 1; i != 0; --i)
            divideStep(other, i);
        return *this;
    }

    template <dirichlet_type Dir> [[nodiscard]] friend Dirichlet operator/(Dirichlet left, const Dir &right)
    {
        left /= right;
        return left;
    }

    /// Multiplication by a scalar.
    Dirichlet &operator*=(T value)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) *= value;
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) *= value;
        return *this;
    }

    [[nodiscard]] friend Dirichlet operator*(Dirichlet S, T value)
    {
        S *= value;
        return S;
    }

    [[nodiscard]] friend Dirichlet operator*(T value, Dirichlet S)
    {
        S *= value;
        return S;
    }

    /// Division by a scalar.
    Dirichlet &operator/=(T value)
    {
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) /= value;
        for (uint32_t i = 1; i < _down.size(); ++i)
            down(i) /= value;
        return *this;
    }

    [[nodiscard]] friend Dirichlet operator/(Dirichlet left, T value)
    {
        left /= value;
        return left;
    }

    /// Dirichlet inverse.
    [[nodiscard]] Dirichlet operator~() const { return inverse(); }

    template <typename U> bool operator==(const Dirichlet<U> &other) const
    {
        return _n == other.n() && _up == other.up() && _down == other.down();
    }

    /// Performs multiplication of S by (f(0) + f(1) * p^-s + ...).
    template <typename Fun> Dirichlet multiplyLocal(this Dirichlet self, size_t p, Fun f)
    {
        self.down() = mapv(std::execution::par, range(0_u32, (uint32_t)self.down().size() - 1),
                           [&](uint32_t i) { return i == 0 ? T(0) : self.localProductValue(p, f, self.quotient(i)); });
        // Sieve for the up values.
        T const c = f(0);
        adjacentDifferenceInPlace(self.up());
        T a{};
        for (size_t i = self.up().size() - 1; i != 0; --i)
        {
            a = self.up(i);
            self.up(i) *= c;
            int e = 1;
            for (size_t j = p; i * j < self.up().size(); ++e, j *= p)
                self.up(i * j) += a * f(e);
        }
        partialSumInPlace(self.up());
        return self;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Dirichlet &S)
    {
        return o << "{\n  n: " << S._n << "\n  up: " << S._up << "\n  down: " << S._down << "\n}";
    }
};

namespace dirichlet
{
/// Computes `∑ (ab ≤ n), f(a) * g(b)` in O(√n).
/// Equivalent to `∑ (k ≤ n), (f * g)(k)`.
/// Equivalent to `∑ (k ≤ n), f(k) * G(n/k) = ∑ (k ≤ n), g(k) * F(n/k)`.
/// Requirements:
/// * `F` and `G` should be the summatory functions of `f` and `g`.
/// * Need to be able to evaluate `f(k)`, `g(k)` for `k ≤ √n` and `F(k)` and `G(k)` for `k ≥ √n`.
template <typename Fun1, typename SummatoryFun1, typename Fun2, typename SummatoryFun2, integral2 T>
auto productValue(Fun1 f, SummatoryFun1 F, Fun2 g, SummatoryFun2 G, T n)
{
    return SpecialDirichlet{std::move(f), std::move(F)}.productValue(SpecialDirichlet{std::move(g), std::move(G)}, n);
}

/// Computes `∑ (ab ≤ n), f(a) * g(b)` in O(√n), given only their summatory functions. Convenience function.
template <typename SummatoryFun1, typename SummatoryFun2, integral2 T>
auto productValue(SummatoryFun1 F, SummatoryFun2 G, T n)
{
    return SpecialDirichlet{[&](auto &&k) { return F(k) - F(k - 1); }, F}.productValue(
        SpecialDirichlet{[&](auto &&k) { return G(k) - G(k - 1); }, G}, n);
}

/// Computes `∑ (p^e * b ≤ n), f(e) * G(b)`.
template <typename Fun, typename SummatoryFun, integral2 T> auto localProductValue(size_t p, Fun f, SummatoryFun G, T n)
{
    using Tp = std::remove_cvref_t<std::invoke_result_t<SummatoryFun, size_t>>;
    assert(p >= 2);
    Tp res = 0;
    for (int e = 0; n != 0; ++e, n /= p)
    {
        auto const c = f(e);
        if (c != 0)
            res += c * G(n);
    }
    return res;
}

/// Computes `∑ (p1^e1 * p2^e2 * b ≤ n), f1(e1) * f2(e2) * G(b)`.
template <typename Fun1, typename Fun2, typename SummatoryFun, integral2 T>
auto localProductValue(size_t p1, Fun1 f1, size_t p2, Fun2 f2, SummatoryFun G, T n)
{
    return localProductValue(p1, f1, [&](size_t k) { return localProductValue(p2, f2, G, k); }, n);
}

/// Computes `∑ (k powerful), f(k) * G(n / k)` in O(√n) evaluations of G.
/// f is given as a multiplicative function: (p, e) ↦ f(p^e).
template <typename Fun, typename SummatoryFun, integral2 T, std::ranges::range Range>
auto powerfulProductValue(Fun f, SummatoryFun G, T n, Range &&primes)
{
    using H = half_integer_t<T>;
    using Tp = std::common_type_t<std::remove_cvref_t<std::invoke_result_t<Fun, T, int>>,
                                  std::remove_cvref_t<std::invoke_result_t<SummatoryFun, T>>>;
    return it::tree(
               std::tuple{T(1), 0UZ, Tp(1)},
               [&](auto &&t, auto rec) {
                   auto &&[k, i, acc] = t;
                   T const hq = fastDiv(n, k);
                   for (size_t j = i; j < primes.size(); ++j)
                   {
                       H const p = primes[j];
                       T const pp = T(p) * p;
                       if (pp > hq)
                           break;
                       int e = 2;
                       for (T q = 1; mulLeq(pp, q, hq); q *= p, ++e)
                       {
                           Tp const value = f(p, e);
                           if (value != 0)
                               rec({k * pp * q, j + 1, acc * value});
                       }
                   }
               },
               [&](auto &&t) {
                   auto &&[k, i, acc] = t;
                   return mulLeq(k, T(primes[i]) * primes[i], n);
               })
        .map([&](auto &&t) {
            auto &&[k, i, acc] = t;
            return acc * G(fastDiv(n, k));
        })
        .sum();
}

/// Computes `∑ (k powerful), f(k) * G(n / k)` in O(√n) evaluations of G.
/// f is given as a multiplicative function: (p, e) ↦ f(p^e).
template <typename Fun, typename SummatoryFun, integral2 T> auto powerfulProductValue(Fun f, SummatoryFun G, T n)
{
    return powerfulProductValue(f, G, n, primeRange<half_integer_t<T>>(isqrt(n)));
}

/// Computes `∑ (k r-powerful), f(k) * G(n / k)` in O(n^(1/r)) evaluations of G.
/// f is given as a multiplicative function: (p, e) ↦ f(p^e).
template <typename Fun, typename SummatoryFun, integral2 T, std::ranges::range Range>
auto powerfulProductValue(Fun f, SummatoryFun G, T n, int r, Range &&primes)
{
    using H = half_integer_t<T>;
    using Tp = std::common_type_t<std::remove_cvref_t<std::invoke_result_t<Fun, T, int>>,
                                  std::remove_cvref_t<std::invoke_result_t<SummatoryFun, T>>>;
    return it::tree(
               std::tuple{T(1), 0UZ, Tp(1)},
               [&](auto &&t, auto rec) {
                   auto &&[k, i, acc] = t;
                   T const hq = fastDiv(n, k);
                   for (size_t j = i; j < primes.size(); ++j)
                   {
                       H const p = primes[j];
                       T const pp = pow(T(p), r);
                       if (pp > hq)
                           break;
                       int e = r;
                       for (T q = 1; mulLeq(pp, q, hq); q *= p, ++e)
                       {
                           Tp const value = f(p, e);
                           if (value != 0)
                               rec({k * pp * q, j + 1, acc * value});
                       }
                   }
               },
               [&](auto &&t) {
                   auto &&[k, i, acc] = t;
                   return mulLeq(k, pow(T(primes[i]), r), n);
               })
        .map([&](auto &&t) {
            auto &&[k, i, acc] = t;
            return acc * G(fastDiv(n, k));
        })
        .sum();
}

/// Computes `∑ (k r-powerful), f(k) * G(n / k)` in O(n^(1/r)) evaluations of G.
/// f is given as a multiplicative function: (p, e) ↦ f(p^e).
template <typename Fun, typename SummatoryFun, integral2 T> auto powerfulProductValue(Fun f, SummatoryFun G, T n, int r)
{
    return powerfulProductValue(f, G, n, r, primeRange<half_integer_t<T>>(inth_root(n, r)));
}

template <typename Fun, integral2 T> inline auto productZeta_2s(Fun &&F, T n)
{
    T const c = inth_root(n, 3);
    return sumMaybeParallel(1, c,
                            [&](T k) { return F(fastDiv(n, k * k)) + isqrt(fastDiv(n, k)) * (F(k) - F(k - 1)); }) -
           c * F(c);
}

template <typename Fun, integral2 T> inline auto productZeta_multiple(int a, Fun &&F, T n)
{
    T const c = inth_root(n, a + 1);
    return sumMaybeParallel(
               1, c, [&](T k) { return F(fastDiv(n, pow(k, a))) + inth_root(fastDiv(n, k), a) * (F(k) - F(k - 1)); }) -
           c * F(c);
}

template <typename Fun, integral2 T> inline auto quotientZeta_multiple(int a, Fun &&F, T n)
{
    assert(a > 1);
    T const c = inth_root(n, a);
    auto const mu = mobiusSieve(c);
    return sumMaybeParallel(1, c, [&](T k) { return mu[k] * F(fastDiv(n, pow(k, a))); });
}

// O(1) Dirichlet series come in two flavors, one is a SpecialDirichlet (recommended) and the other is a Dirichlet.

/// 1, the multiplicative identity. f(n) = [n = 1]. Motive = 0.
template <typename T = int64_t> constexpr auto unit()
{
    return SpecialDirichlet{[](size_t k) -> T { return k == 1; }, [](size_t) -> T { return 1; }};
}

/// 1, the multiplicative identity. f(n) = [n = 1]. Motive = 0.
template <typename T = int64_t> Dirichlet<T> unit(size_t n)
{
    return {n, [](size_t) -> T { return 1; }};
}

/// ζ(s). f(n) = 1. Motive = [1].
template <typename T = int64_t> constexpr auto zeta()
{
    return SpecialDirichlet{[](size_t) -> T { return 1; }, [](size_t k) -> T { return k; }};
}

/// ζ(s). f(n) = 1. Motive = [1].
template <typename T = int64_t> Dirichlet<T> zeta(size_t n) { return {n, std::identity{}}; }

/// ζ(s - 1). f(n) = n. Motive = [p].
template <typename T = int64_t> constexpr auto id()
{
    return SpecialDirichlet{[](size_t k) -> T { return k; }, [](size_t k) -> T { return sumId<T>(k); }};
}

/// ζ(s - 1). f(n) = n. Motive = [p].
template <typename T = int64_t> Dirichlet<T> id(size_t n)
{
    return {n, [](size_t k) -> T { return sumId<T>(k); }};
}

/// ζ(s - 2). f(n) = n^2. If `T = ZMod<M>`, currently requires that neither 2 or 3 divide `M`. Motive = [p^2].
template <typename T = int64_t> constexpr auto id2()
{
    return SpecialDirichlet{[](size_t k) -> T { return T(k) * T(k); }, [](size_t k) -> T { return sumSquares<T>(k); }};
}

/// ζ(s - 2). f(n) = n^2. If `T = ZMod<M>`, currently requires that neither 2 or 3 divide `M`. Motive = [p^2].
template <typename T = int64_t> Dirichlet<T> id2(size_t n)
{
    // Fancy stuff to avoid overflow or ZMod divisions.
    return {n, [](size_t k) -> T { return sumSquares<T>(k); }};
}

/// ζ(2s). f(n) = [n is square]. Motive = [1] + [-1].
template <typename T = int64_t> constexpr auto zeta_2s()
{
    return SpecialDirichlet{[](size_t k) -> T { return isSquare(k); }, [](size_t k) -> T { return isqrt(k); }};
}

/// ζ(2s). f(n) = [n is square]. Motive = [1] + [-1].
template <typename T = int64_t> Dirichlet<T> zeta_2s(size_t n)
{
    return {n, [](size_t k) -> T { return isqrt(k); }};
}

// TODO: Add SpecialDirichlets for these at some point.

/// ζ(s - 3). f(n) = n^3. Motive = [p^3].
template <typename T = int64_t> Dirichlet<T> id3(size_t n)
{
    return {n, [](size_t k) -> T {
                T const s = sumId<T>(k);
                return s * s;
            }};
}

/// ζ(as). f(n) = [n is a perfect ath power]. Motive = sum([r] for r in ath roots of unity)
template <typename T = int64_t> Dirichlet<T> zeta_multiple(int a, size_t n)
{
    return {n, [&](size_t k) -> T { return inth_root(k, a); }};
}

/// χ_(-3)(s). f(n) = (-3|n). 1 if 1 mod 3, -1 if 2 mod 3, 0 otherwise. L-function of the hexagonal lattice and the
/// Eisenstein integers.
template <typename T = int64_t> Dirichlet<T> chi3(size_t n)
{
    return {n, [&](size_t k) -> T { return k % 3 == 1; }};
}

/// χ_(-4)(s). f(n) = (-4|n). 1 if 1 mod 4, -1 if 3 mod 4, 0 otherwise. Also known as the Dirichlet beta function.
/// Motive = χ_(-4).
template <typename T = int64_t> Dirichlet<T> chi4(size_t n)
{
    return {n, [&](size_t k) -> T { return k % 4 == 1 || k % 4 == 2; }};
}

/// χ_5(s). f(n) = (5|n).
template <typename T = int64_t> Dirichlet<T> chi5(size_t n)
{
    return {n, [&](size_t k) -> T { return std::array{0, 1, 0, -1, 0}[k % 5]; }};
}

/// 1 / ζ(2s). f(n) = [n is square] * μ(√n). O(n^(1/2)). Motive = -[1] - [-1].
template <typename T = int64_t> Dirichlet<T> inv_zeta_2s(size_t n)
{
    auto const mu = mobiusSieve(isqrt(n));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](size_t k) -> T { return mertens[isqrt(k)]; }};
}

/// 1 / ζ(as). f(n) = [n is a perfect rth power] * μ(n^(1/a)). O(n^(1/a)). Motive = -∑([r] for r in ath roots of unity).
template <typename T = int64_t> Dirichlet<T> inv_zeta_multiple(int a, size_t n)
{
    auto const mu = mobiusSieve(inth_root(n, a));
    auto const mertens = partialSum(mu, T{});
    return {n, [&](size_t k) -> T { return mertens[inth_root(k, a)]; }};
}

/// ζ(as - b). f(n) = [n = k^a] * k^b. Requires `a > 1`. Motive = ∑([r*p^(b/a)] for r in ath roots of unity).
template <typename T = int64_t> Dirichlet<T> zeta_linear(int a, int b, size_t n)
{
    assert(a > 1);
    size_t const s = inth_root(n, a);
    auto sieve = range(0, s, [&](size_t k) { return std::pow(T(k), b); });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](size_t k) -> T { return sieve[inth_root(k, a)]; }};
}

/// 1 / ζ(as - b). f(n) = [n = k^a] * μ(k) * k^b. Requires `a > 1`. Motive = -∑([r*p^(b/a)] for r in ath roots of
/// unity).
template <typename T = int64_t> Dirichlet<T> inv_zeta_linear(int a, int b, size_t n)
{
    assert(a > 1);
    size_t const s = inth_root(n, a);
    auto const mu = mobiusSieve(s);
    auto sieve = range(0, s, [&](size_t k) { return std::pow(T(k), b) * mu[k]; });
    sieve[0] = 0;
    partialSumInPlace(sieve);
    return {n, [&](size_t k) -> T { return sieve[inth_root(k, a)]; }};
}

/// ζ(s) / ζ(2s). f(n) = |μ(n)| = [n is squarefree]. O(n^(3/5)). Motive = -[-1].
template <typename T = int64_t> Dirichlet<T> squarefree(size_t n, double alpha = 1.25)
{
    size_t const s =
        std::min(Dirichlet<T>::pivotMax, std::max(Dirichlet<T>::defaultPivot(n), (size_t)(alpha * std::pow(n, 0.6))));
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
template <typename T = int64_t> Dirichlet<T> tau(size_t n, double alpha = 0.08)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivotMax, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    auto const precomputed = partialSum(divisorCountSieve(s), T{});
    Dirichlet<T> S{n};
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
template <typename T = int64_t> Dirichlet<T> sigma(size_t n, double alpha = 0.08)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivotMax, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    return id<T>(n).multiply(zeta<T>(), partialSum(divisorSumSieve(s), T{}));
}

/// ζ(s)ζ(s - 2). f(n) = sum of squares of divisors of n. O(n^(2/3)). Motive = [p^2] + [1].
template <typename T = int64_t> Dirichlet<T> sigma2(size_t n) { return id2<T>(n).multiply(zeta<T>()); }

/// ζ(s)ζ(s - 3). f(n) = sum of cubes of divisors of n. O(n^(2/3)). Motive = [p^3] + [1].
template <typename T = int64_t> Dirichlet<T> sigma3(size_t n) { return id3<T>(n).multiply(zeta<T>()); }

/// 1 / ζ(s). f(n) = μ(n). F(n) is the Mertens function. O(n^(2/3)). Motive = -[1].
template <typename T = int64_t> Dirichlet<T> mobius(size_t n, double alpha = 0.45)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivotMax, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    return unit<T>(n).divide(zeta<T>(), partialSum(mobiusSieve(s), T{}));
}

/// ζ(s - 1) / ζ(s). f(n) = φ(n). O(n^(2/3)). Motive = [p] - [1].
template <typename T = int64_t> Dirichlet<T> totient(size_t n, double alpha = 0.35)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivotMax, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    return id<T>(n).divide(zeta<T>(), partialSum(totientSieve(s), T{}));
}

/// ζ(2s) / ζ(s). f(n) = (-1)^(number of primes dividing n). O(n^(2/3)). Motive = [-1].
template <typename T = int64_t> Dirichlet<T> liouville(size_t n, double alpha = 0.35)
{
    size_t const s = std::max(Dirichlet<T>::defaultPivot(n),
                              std::min(Dirichlet<T>::pivotMax, (size_t)(alpha * std::pow(n, 2.0 / 3))));
    return zeta_2s<T>(n).divide(zeta<T>(), partialSum(liouvilleSieve(s), T{}));
}
} // namespace dirichlet
} // namespace euler
