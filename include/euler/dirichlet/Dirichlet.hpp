#pragma once

#include "../it/tree.hpp"
#include "libdivide.h"
#include "sieve_ops.hpp"

namespace euler
{
template <typename T = i64> class Dirichlet;
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
    Fun f_;
    SFun F_;

  public:
    constexpr SpecialDirichlet(Fun f, SFun F) : f_(std::move(f)), F_(std::move(F)) {}

    [[nodiscard]] constexpr auto value(size_t k) const noexcept { return f_(k); }
    [[nodiscard]] constexpr auto sum(size_t k) const noexcept { return F_(k); }
    [[nodiscard]] constexpr auto up(size_t k) const noexcept { return F_(k); }

    [[nodiscard]] constexpr auto operator[](size_t k) const noexcept { return F_(k); }
    [[nodiscard]] constexpr auto operator()(size_t k) const noexcept { return F_(k); }

    template <typename Fun2, typename SFun2, integral2 T>
    [[nodiscard]] auto productValue(const SpecialDirichlet<Fun2, SFun2> &other, T n) const
    {
        using H = half_integer_t<T>;
        using Tp = std::common_type_t<std::invoke_result_t<SFun, T>, std::invoke_result_t<SFun2, T>>;
        H const s = (H)isqrt(n);
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
    size_t n_ = 0;
    std::vector<T> up_;
    std::vector<T> down_;
    std::vector<size_t> quots_;

    Dirichlet(size_t n, std::vector<T> up, std::vector<T> down, std::vector<size_t> quots)
        : n_(n), up_(std::move(up)), down_(std::move(down)), quots_(std::move(quots))
    {
    }

    template <dirichlet_type Dir> T sumTerm1(const Dir &other, size_t i, u32 j) const
    {
        if constexpr (requires { other.down(i * j); })
            return value(j) * other.down(i * j) + down_[i * j] * other.value(j);
        return value(j) * other.sum(quotient(i * j)) + down_[i * j] * other.value(j);
    }

    template <dirichlet_type Dir> T sumTerm2(const Dir &other, const libdivide::divider<size_t> &i, u32 j) const
    {
        size_t const k = quotient(j) / i;
        return value(j) * other.up(k) + up_[k] * other.value(j);
    }

    /// One step of the divison algorithm by ζ(s). Internal use only.
    template <dirichlet_type Dir> void divideStep(const Dir &other, u32 i)
    {
        libdivide::divider<size_t> const fast_i(i);
        size_t const k = quotient(i);
        u32 const s = isqrt(k);
        u32 const mid = (down_.size() - 1) / i;
        down_[i] -= (up(1) - up(0)) * T(other.sum(k));
        for (u32 j = 2; j <= mid; ++j)
            down_[i] -= sumTerm1(other, i, j);
        for (u32 j = mid + 1; j <= s; ++j)
            down_[i] -= sumTerm2(other, fast_i, j);
        down_[i] += up(s) * T(other.up(s));
        if (other.value(1) != 1)
            down_[i] /= other.value(1);
    }

  public:
    Dirichlet() = default;

    /// Default max memory usage = 24 GB.
    static constexpr size_t DefaultMaxMemoryUsage = 3UZ << 33;
    /// Configurable upper bound on the size of the up vector.
    static inline size_t pivot_max = DefaultMaxMemoryUsage / sizeof(T);
    /// Configurable exponent on the size of the up vector.
    static inline double pivot_exponent = 2.0 / 3;
    /// Configurable coefficient on the size of the up vector. Ideal is 0.2 for multiplication and 0.5 for division.
    static inline double pivot_coefficient = 1;

    /// Gives the default size of the up vector for a given value of `n`.
    static size_t defaultPivot(size_t n)
    {
        size_t const s = pivot_coefficient * std::pow(n / std::max(1.0, log(n)), pivot_exponent);
        // Impose maximum so as to not overwhelm memory.
        size_t const res = std::max(isqrt(n), std::min(pivot_max, s));
        return n / (n / res);
    }

    explicit Dirichlet(size_t n)
        : n_(n), up_(defaultPivot(n_) + 1), down_(n_ / (defaultPivot(n_) + 1) + 1), quots_(isqrt(n_) + 1)
    {
        for (size_t i = 1; i < quots_.size(); ++i)
            quots_[i] = n_ / i;
    }

    /// Initialize with the given summatory function.
    template <typename SummatoryFun> Dirichlet(size_t n, SummatoryFun F) : Dirichlet(n)
    {
        for (size_t k = 1; k < up_.size(); ++k)
            up_[k] = T(F(k));
        for (u32 i = down_.size() - 1; i != 0; --i)
            down_[i] = T(F(quots_[i]));
    }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] size_t n() const noexcept { return n_; }
    /// The index of the transition point between up and down vectors.
    [[nodiscard]] size_t pivot() const noexcept { return up_.size() - 1; }
    /// Gets the value of `n / i` avoiding a division CPU instruction. `i` must be `≤ √n`.
    [[nodiscard]] size_t quotient(u32 i) const noexcept { return quots_[i]; }

    /// The up vector, mutable.
    [[nodiscard]] std::vector<T> &up() noexcept { return up_; }
    /// The up vector.
    [[nodiscard]] const std::vector<T> &up() const noexcept { return up_; }
    /// Element access into the up array, mutable.
    [[nodiscard]] T &up(size_t x) noexcept { return up_[x]; }
    /// Element access into the up array.
    [[nodiscard]] const T &up(size_t x) const noexcept { return up_[x]; }

    /// The down vector, mutable.
    [[nodiscard]] std::vector<T> &down() noexcept { return down_; }
    /// The down vector.
    [[nodiscard]] const std::vector<T> &down() const noexcept { return down_; }
    /// Element access into the down array, mutable.
    [[nodiscard]] T &down(size_t x) noexcept { return down_[x]; }
    /// Element access into the down array.
    [[nodiscard]] const T &down(size_t x) const noexcept { return down_[x]; }

    [[nodiscard]] T &front() noexcept { return up(1); }
    [[nodiscard]] const T &front() const noexcept { return up(1); }
    [[nodiscard]] T &back() noexcept { return down_.size() > 1 ? down(1) : up_.back(); }
    [[nodiscard]] const T &back() const noexcept { return down_.size() > 1 ? down(1) : up_.back(); }

    /// Returns a vector of function values from 1 to the pivot. This is the vector of adjacent differences of `up`.
    [[nodiscard]] std::vector<T> values() const { return adjacentDifference(up_); }
    /// Returns a single value. Requires `k ≤ pivot()`.
    [[nodiscard]] T value(size_t k) const noexcept { return up(k) - up(k - 1); }

    [[nodiscard]] T &operator[](size_t k) noexcept { return k < up_.size() ? up_[k] : down_[n_ / k]; }
    [[nodiscard]] const T &operator[](size_t k) const noexcept { return k < up_.size() ? up_[k] : down_[n_ / k]; }
    [[nodiscard]] T &operator()(size_t k) noexcept { return (*this)[k]; }
    [[nodiscard]] const T &operator()(size_t k) const noexcept { return (*this)[k]; }

    /// Returns the prefix sum of this series at the given index.
    [[nodiscard]] const T &sum(size_t k) const noexcept { return (*this)[k]; }

    /// Takes a partial sum in place. Useful function to call when building this array from individual values.
    Dirichlet &accumulate()
    {
        for (size_t i = 1; i < up_.size(); ++i)
            up(i) += up(i - 1);
        down_.back() += up_.back();
        for (u32 i = down_.size() - 2; i != 0; --i)
            down(i) += down(i + 1);
        return *this;
    }

    /// Applies a function to each entry of this array.
    template <std::invocable<T> Fun> Dirichlet &mapInPlace(Fun f)
    {
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) = f(up(k));
        for (u32 i = down_.size() - 1; i != 0; --i)
            down(i) = f(down(i));
        return *this;
    }

    /// Enumerates keys of this floors array in ascending order. Breaks if `f` returns `it::result_break`.
    template <std::invocable<size_t> Fun> it::result_t ascending(Fun f) const
    {
        for (size_t k = 1; k < up_.size(); ++k)
            if (!it::callbackResult(f, k))
                return it::result_break;
        for (u32 i = down_.size() - 1; i != 0; --i)
            if (!it::callbackResult(f, quotient(i)))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (keys, mutable value) pairs of this floors array in ascending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> it::result_t ascendingMut(Fun f)
    {
        for (size_t k = 1; k < up_.size(); ++k)
            if (!it::callbackResult(f, k, up(k)))
                return it::result_break;
        for (u32 i = down_.size() - 1; i != 0; --i)
            if (!it::callbackResult(f, quotient(i), down(i)))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates keys of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t> Fun> it::result_t descending(Fun f) const
    {
        for (u32 i = 1; i < down_.size(); ++i)
            if (!it::callbackResult(f, quotient(i)))
                return it::result_break;
        for (size_t k = up_.size() - 1; k != 0; --k)
            if (!it::callbackResult(f, k))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (key, mutable value) pairs of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> it::result_t descendingMut(Fun f)
    {
        for (u32 i = 1; i < down_.size(); ++i)
            if (!it::callbackResult(f, quotient(i), down(i)))
                return it::result_break;
        for (size_t k = up_.size() - 1; k != 0; --k)
            if (!it::callbackResult(f, k, up(k)))
                return it::result_break;
        return it::result_continue;
    }

    /// Compute a single summatory value of `this * this`.
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <size_t ParThreshold = 8192> [[nodiscard]] T squareValue(size_t k) const
    {
        u32 const s = isqrt(k);
        size_t const i = n_ / k;
        u32 const u = (down_.size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return 2 * (sumMaybeParallel<ParThreshold>(1, u, [&](u32 j) -> T { return value(j) * down(i * j); }) +
                    sumMaybeParallel<ParThreshold>(u + 1, s,
                                                   [&](u32 j) -> T { return value(j) * up(quotient(j) / fasti); })) -
               up(s) * up(s);
    }

    /// Compute a single summatory value of (this * other).
    /// Precondition: `k` must be of the form `⌊n / i⌋` for some `i`.
    template <size_t ParThreshold = 8192, dirichlet_type Dir>
    [[nodiscard]] T productValue(const Dir &other, size_t k) const
    {
        if constexpr (std::is_same_v<std::decay_t<Dir>, Dirichlet>)
            if (this == &other)
                return squareValue(k);
        u32 const s = isqrt(k);
        size_t const i = n_ / k;
        u32 const u = (down().size() - 1) / i;
        libdivide::divider<size_t> const fasti(i);
        return sumMaybeParallel<ParThreshold>(1, u, [&](u32 j) -> T { return sumTerm1(other, i, j); }) +
               sumMaybeParallel<ParThreshold>(u + 1, s, [&](u32 j) -> T { return sumTerm2(other, fasti, j); }) -
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

    /// Dirichlet inverse.
    [[nodiscard]] Dirichlet inverse() const
    {
        Dirichlet res{n_, [](auto &&) { return T(1); }};
        res /= *this;
        return res;
    }

    /// Returns this raised to a power.
    template <integral2 E> [[nodiscard]] Dirichlet pow(this Dirichlet self, E exponent)
    {
        if (exponent == 0)
            return Dirichlet{self.n(), [](auto &&) { return T(1); }};
        if (exponent < 0)
        {
            self = ~self;
            exponent = -exponent;
        }
        for (; !(exponent & 1); exponent >>= 1)
            self.squareInPlace();
        if (exponent == 1)
            return self;
        Dirichlet x = self;
        for (exponent >>= 1; exponent != 0; exponent >>= 1)
        {
            self.squareInPlace();
            if (exponent & 1)
                x *= self;
        }
        return x;
    }

    /// Addition.
    template <typename U> Dirichlet &operator+=(const Dirichlet<U> &other)
    {
        assert(n_ == other.n());
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) += other.up(k);
        for (u32 i = 1; i < down_.size(); ++i)
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
        assert(n_ == other.n());
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) -= other.up(k);
        for (u32 i = 1; i < down_.size(); ++i)
            down(i) -= other.down(i);
        return *this;
    }

    template <typename U> [[nodiscard]] friend Dirichlet operator-(Dirichlet left, const Dirichlet<U> &right)
    {
        left -= right;
        return left;
    }

    /// Squares this Dirichlet series in place.
    /// Precondition: `precomp` must be at least as large as the up vector, or be empty.
    template <std::ranges::range Range = std::ranges::empty_view<T>> Dirichlet &squareInPlace(Range &&precomp = {})
    {
        u32 const max_i = precomp.empty() ? down_.size() - 1 : n_ / precomp.size();
        down_ = mapv(std::execution::par, std::views::iota(0U, (u32)down_.size()),
                     [&](u32 i) { return i == 0 || i > max_i ? T(0) : squareValue<0>(quotient(i)); });

        if (precomp.empty())
        {
            // Sieve for the up values.
            adjacentDifferenceInPlace(up_);
            dirichlet::convolveInPlace(std::execution::par, up_, up_);
            partialSumInPlace(std::execution::par, up_);
        }
        else
        {
            for (u32 i = max_i + 1; i < down_.size(); ++i)
                down(i) = precomp[quotient(i)];
            std::copy_n(std::execution::par, precomp.begin(), up_.size(), up_.begin());
        }
        return *this;
    }

    /// Precondition: `precomp` must be at least as large as the up vector, or be empty.
    template <std::ranges::range Range = std::ranges::empty_view<T>>
    [[nodiscard]] Dirichlet square(this Dirichlet self, Range &&precomp = {})
    {
        self.squareInPlace(std::forward<Range>(precomp));
        return self;
    }

    /// Multiplication in place with a precomputed sieve.
    /// Precondition: `precomp` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    Dirichlet &multiplyInPlace(const Dir &other, Range &&precomp = {})
    {
        if (precomp.empty())
            return *this *= other;
        u32 const max_i = n_ / precomp.size();
        down_ = mapv(std::execution::par, std::views::iota(0U, (u32)down_.size()),
                     [&](u32 i) { return i == 0 || i > max_i ? T(0) : productValue<0>(other, quotient(i)); });
        for (u32 i = max_i + 1; i < down_.size(); ++i)
            down_[i] = precomp[quots_[i]];
        std::copy_n(std::execution::par, precomp.begin(), up_.size(), up_.begin());
        return *this;
    }

    /// Multiplication with a precomputed sieve.
    /// Precondition: `precomp` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    [[nodiscard]] Dirichlet multiply(this Dirichlet self, const Dir &other, Range &&precomp = {})
    {
        self.multiplyInPlace(other, std::forward<Range>(precomp));
        return self;
    }

    /// Multiplication.
    template <dirichlet_type Dir> Dirichlet &operator*=(const Dir &other)
    {
        if constexpr (std::is_same_v<std::decay_t<Dir>, Dirichlet>)
            if (this == &other)
                return squareInPlace();
        for (u32 start = 1; start < down_.size(); start <<= 1)
            tbb::parallel_for(start, std::min((u32)down_.size(), start << 1),
                              [&](u32 i) { down_[i] = productValue<0>(other, quots_[i]); });
        adjacentDifferenceInPlace(up_);
        dirichlet::convolveInPlace(std::execution::par, up_, [&](size_t i) { return other.value(i); });
        partialSumInPlace(std::execution::par, up_);
        return *this;
    }

    /// Multiplication.
    template <dirichlet_type Dir> [[nodiscard]] Dirichlet operator*(this Dirichlet self, const Dir &other)
    {
        self *= other;
        return self;
    }

    template <typename Fun, typename SFun>
    [[nodiscard]] friend Dirichlet operator*(const SpecialDirichlet<Fun, SFun> &left, const Dirichlet &right)
    {
        return right * left;
    }

    /// Division in place with a precomputed sieve.
    /// Precondition: `precomp` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    Dirichlet &divideInPlace(const Dir &other, Range &&precomp = {})
    {
        if (precomp.empty())
            return *this /= other;
        std::copy_n(std::execution::par, precomp.begin(), up_.size(), up_.begin());
        u32 max_i = n_ / precomp.size();
        for (u32 i = max_i + 1; i < down_.size(); ++i)
            down_[i] = precomp[quots_[i]];
        for (; max_i != 0; max_i >>= 1)
            tbb::parallel_for((max_i >> 1) + 1, max_i + 1, [&](u32 i) { divideStep(other, i); });
        return *this;
    }

    /// Division with a precomputed sieve.
    /// Precondition: `precomp` must be at least as large as the up vector, or be empty.
    template <dirichlet_type Dir, std::ranges::range Range = std::ranges::empty_view<T>>
    [[nodiscard]] Dirichlet divide(this Dirichlet self, const Dir &other, Range &&precomp = {})
    {
        self.divideInPlace(other, std::forward<Range>(precomp));
        return self;
    }

    /// Division.
    template <dirichlet_type Dir> Dirichlet &operator/=(const Dir &other)
    {
        auto const c = other.value(1);
        // Sieve for the up values.
        adjacentDifferenceInPlace(up_);
        dirichlet::invConvolve(std::execution::par, up_, [&](size_t i) { return other.value(i); });
        partialSumInPlace(std::execution::par, up_);
        for (u32 max_i = down_.size() - 1; max_i != 0; max_i >>= 1)
            tbb::parallel_for((max_i >> 1) + 1, max_i + 1, [&](u32 i) { divideStep(other, i); });
        return *this;
    }

    template <dirichlet_type Dir> [[nodiscard]] Dirichlet operator/(this Dirichlet self, const Dir &other)
    {
        self /= other;
        return self;
    }

    /// Multiplication by a scalar.
    Dirichlet &operator*=(T value)
    {
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) *= value;
        for (u32 i = 1; i < down_.size(); ++i)
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
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) /= value;
        for (u32 i = 1; i < down_.size(); ++i)
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
        return n_ == other.n() && up_ == other.up() && down_ == other.down();
    }

    /// Performs multiplication of S by (f(0) + f(1) * p^-s + ...).
    template <typename Fun> Dirichlet multiplyLocal(this Dirichlet self, size_t p, Fun f)
    {
        self.down() = mapv(std::execution::par, std::views::iota(0U, (u32)self.down().size()),
                           [&](u32 i) { return i == 0 ? T(0) : self.localProductValue(p, f, self.quotient(i)); });
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
        partialSumInPlace(std::execution::par, self.up());
        return self;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Dirichlet &S)
    {
        return o << "{\n  n: " << S.n_ << "\n  up: " << S.up_ << "\n  down: " << S.down_ << "\n}";
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
} // namespace dirichlet
} // namespace euler
