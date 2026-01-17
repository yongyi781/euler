#pragma once

#include "decls.hpp"
#include "it/primes.hpp"
#include "libdivide.h"
#include "prime.hpp"
#include <ranges>

namespace euler
{
/// An efficient sparse array whose keys are floors ⌊N/i⌋ for i >= 1.
/// Usage: `floors_array(n)` or `floors_array(n, s)`.
template <typename T = int64_t> class floors_array
{
    size_t n_ = 0;
    std::vector<T> up_;
    std::vector<T> down_;
    std::vector<size_t> quots_;

  public:
    /// s should be 1 less than the size of the small values array.
    constexpr floors_array() = default;
    constexpr floors_array(size_t n, size_t pivot)
        : n_(n), up_(pivot + 1), down_(n_ / (pivot + 1) + 1), quots_(isqrt(n_) + 1)
    {
        for (size_t i = 1; i < quots_.size(); ++i)
            quots_[i] = n_ / i;
    }
    constexpr explicit floors_array(size_t n) : floors_array(n, isqrt(n)) {}

    constexpr T &operator[](size_t i) { return i < up_.size() ? up_[i] : down_[n_ / i]; }
    constexpr const T &operator[](size_t i) const { return i < up_.size() ? up_[i] : down_[n_ / i]; }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] constexpr size_t n() const { return n_; }
    /// The transition point between up and down. What was passed as the s parameter during construction.
    [[nodiscard]] constexpr size_t pivot() const { return up_.size() - 1; }
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

    [[nodiscard]] constexpr T &front() { return up_[1]; }
    [[nodiscard]] constexpr const T &front() const { return up_[1]; }
    [[nodiscard]] constexpr T &back() { return down_[1]; }
    [[nodiscard]] constexpr const T &back() const { return down_[1]; }

    /// Addition.
    template <typename U> floors_array &operator+=(const floors_array<U> &other)
    {
        assert(n_ == other.n());
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) += other.up(k);
        for (u32 i = 1; i < down_.size(); ++i)
            down(i) += other.down(i);
        return *this;
    }

    template <typename U> [[nodiscard]] friend floors_array operator+(floors_array left, const floors_array<U> &right)
    {
        left += right;
        return left;
    }

    /// Subtraction.
    template <typename U> floors_array &operator-=(const floors_array<U> &other)
    {
        assert(n_ == other.n());
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) -= other.up(k);
        for (u32 i = 1; i < down_.size(); ++i)
            down(i) -= other.down(i);
        return *this;
    }

    template <typename U> [[nodiscard]] friend floors_array operator-(floors_array left, const floors_array<U> &right)
    {
        left -= right;
        return left;
    }

    /// Division by a scalar.
    floors_array &operator/=(T value)
    {
        for (size_t k = 1; k < up_.size(); ++k)
            up(k) /= value;
        for (u32 i = 1; i < down_.size(); ++i)
            down(i) /= value;
        return *this;
    }

    [[nodiscard]] friend floors_array operator/(floors_array left, T value)
    {
        left /= value;
        return left;
    }

    /// Enumerates keys of this floors array in ascending order. Breaks if `f` returns `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t ascending(Fun f) const
    {
        for (size_t i = 1; i < up_.size(); ++i)
            if (!it::callbackResult(f, i))
                return it::result_break;
        for (size_t i = down_.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, quots_[i]))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (keys, mutable value) pairs of this floors array in ascending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t ascendingMut(Fun f)
    {
        for (size_t i = 1; i < up_.size(); ++i)
            if (!it::callbackResult(f, i, up_[i]))
                return it::result_break;
        for (size_t i = down_.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, quots_[i], down_[i]))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates keys of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t> Fun> constexpr it::result_t descending(Fun f) const
    {
        for (size_t i = 1; i < down_.size(); ++i)
            if (!it::callbackResult(f, quots_[i]))
                return it::result_break;
        for (size_t i = up_.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, i))
                return it::result_break;
        return it::result_continue;
    }

    /// Enumerates (key, mutable value) pairs of this floors array in descending order. Breaks if `f` returns
    /// `it::result_break`.
    template <std::invocable<size_t, T &> Fun> constexpr it::result_t descendingMut(Fun f)
    {
        for (size_t i = 1; i < down_.size(); ++i)
            if (!it::callbackResult(f, quots_[i], down_[i]))
                return it::result_break;
        for (size_t i = up_.size() - 1; i > 0; --i)
            if (!it::callbackResult(f, i, up_[i]))
                return it::result_break;
        return it::result_continue;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const floors_array &S)
    {
        return o << "{\n  n: " << S.n_ << "\n  up: " << S.up_ << "\n  down: " << S.down_ << "\n}";
    }
};

/// Returns a floors array of values `(1 ≤ p ≤ k, p prime) * f(p)` for `k` of the form `⌊limit / i⌋`. Here, `f` must be
/// a completely multiplicative function and `F` must be the summatory function of `f`.
template <typename Fun, typename SummatoryFun> auto primeSumTable(size_t N, Fun f, SummatoryFun F)
{
    using T = std::invoke_result_t<SummatoryFun, size_t>;
    floors_array<T> res(N);
    u32 const s = isqrt(N);
    res.ascendingMut([&](size_t i, T &value) { value = F(i) - F(1); });
    bool const use_omp = N >= 50'000'000'000UZ;
    it::primes(2, s)([&](u64 p) {
        T const fp = f(p);
        T const sp_prev = res.up(p - 1);
        size_t const pp = p * p;
        size_t const mid_i = (res.down().size() - 1) / p;
        size_t const max_i = std::min(res.down().size() - 1, N / pp);
        for (size_t i = 1; i <= mid_i; ++i)
            res.down(i) -= fp * (res.down(i * p) - sp_prev);
        libdivide::divider const fastp(p);
        if (use_omp)
#pragma omp parallel for schedule(static) if (max_i > mid_i + 2048)
            for (size_t i = mid_i + 1; i <= max_i; ++i)
                res.down(i) -= fp * (res.up(res.quotient(i) / fastp) - sp_prev);
        else
            for (size_t i = mid_i + 1; i <= max_i; ++i)
                res.down(i) -= fp * (res.up(res.quotient(i) / fastp) - sp_prev);
        if (pp < res.up().size())
            for (size_t k = res.up().size() - 1, q = k / fastp; k >= pp; --q)
            {
                T const val = fp * (res.up(q) - sp_prev);
                size_t const min_k = std::max(pp, q * p);
                for (; k >= min_k; --k)
                    res.up(k) -= val;
            }
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
