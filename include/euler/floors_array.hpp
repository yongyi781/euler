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
    size_t _n = 0;
    std::vector<T> _up;
    std::vector<T> _down;

  public:
    /// s should be 1 less than the size of the small values array.
    constexpr floors_array() = default;
    constexpr floors_array(size_t n, size_t pivot) : _n(n), _up(pivot + 1), _down(n / (pivot + 1) + 1) {}
    constexpr explicit floors_array(size_t n) : floors_array(n, isqrt(n)) {}

    constexpr T &operator[](size_t i) { return i < _up.size() ? _up[i] : _down[_n / i]; }
    constexpr const T &operator[](size_t i) const { return i < _up.size() ? _up[i] : _down[_n / i]; }

    /// The number that this array was designed for, i.e. the top index.
    [[nodiscard]] constexpr size_t n() const { return _n; }
    /// The transition point between up and down. What was passed as the s parameter during construction.
    [[nodiscard]] constexpr size_t pivot() const { return _up.size() - 1; }

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

    [[nodiscard]] constexpr T &front() { return _up[1]; }
    [[nodiscard]] constexpr const T &front() const { return _up[1]; }
    [[nodiscard]] constexpr T &back() { return _down[1]; }
    [[nodiscard]] constexpr const T &back() const { return _down[1]; }

    /// Addition.
    template <typename U> floors_array &operator+=(const floors_array<U> &other)
    {
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) += other.up(k);
        for (uint32_t i = 1; i < _down.size(); ++i)
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
        assert(_n == other.n());
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) -= other.up(k);
        for (uint32_t i = 1; i < _down.size(); ++i)
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
        for (size_t k = 1; k < _up.size(); ++k)
            up(k) /= value;
        for (uint32_t i = 1; i < _down.size(); ++i)
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

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const floors_array &S)
    {
        return o << "{\n  n: " << S._n << "\n  up: " << S._up << "\n  down: " << S._down << "\n}";
    }
};

/// Returns a floors array of values `(1 ≤ p ≤ k, p prime) * f(p)` for `k` of the form `⌊limit / i⌋`. Here, `f` must be
/// a completely multiplicative function and `F` must be the summatory function of `f`.
template <std::invocable<size_t> Fun, std::invocable<size_t> SummatoryFun>
constexpr auto primeSumTable(size_t limit, Fun f, SummatoryFun F)
{
    using FT = std::remove_cvref_t<std::invoke_result_t<SummatoryFun, int64_t>>;
    size_t const r = isqrt(limit);
    floors_array<FT> S(limit);
    auto const invs = range(0UZ, r, [&](size_t i) { return i == 0 ? 0 : limit / i; });
    S.ascendingMut([&](size_t i, FT &value) { value = F(i) - F(1); });
    for (size_t p = 2; p <= r; ++p)
    {
        if (S.up()[p] == S.up()[p - 1])
            continue;
        // p is prime at this point.
        auto const fp = f(p);
        size_t const pp = p * p;
        size_t const hi1 = (S.down().size() - 1) / p;
        size_t const hi = std::min(S.down().size() - 1, limit / pp);
        // Iterate downward in 3 steps. First from top to p * √N.
        for (size_t i = 1; i <= hi1; ++i)
            S.down()[i] -= fp * (S.down()[i * p] - S.up()[p - 1]);
        // Second from p * √N to √N.
        libdivide::divider const fastp(p);
        for (size_t i = hi1 + 1; i <= hi; ++i)
            S.down()[i] -= fp * (S.up()[invs[i] / fastp] - S.up()[p - 1]);
        if (pp <= S.up().size() - 1)
        {
            // Finally from √N to 1.
            for (size_t k = S.up().size() - 1; k >= pp; --k)
                S.up()[k] -= fp * (S.up()[k / fastp] - S.up()[p - 1]);
        }
    }
    return S;
}

/// Returns a floors array of values `(1 ≤ p ≤ k, p prime) * p` for `k` of the form `⌊limit / i⌋`.
template <typename T = int64_t> constexpr floors_array<T> primeSumTable(size_t limit)
{
    return primeSumTable(limit, [](size_t n) -> T { return n; }, [](size_t n) -> T { return sumId<T>(n); });
}

/// Returns a floors array of values `#(1 ≤ p ≤ k, p prime)` for `k` of the form `⌊limit / i⌋`.
template <typename T = int64_t> constexpr floors_array<T> primePiTable(size_t limit)
{
    return primeSumTable(limit, [](size_t) -> T { return 1; }, [](size_t n) -> T { return n; });
}

/// Calculates `(1 ≤ p ≤ limit, p prime) * f(p)`. Here, `f` must be a completely multiplicative function and `F` must be
/// the summatory function of `f`.
template <std::invocable<int64_t> Fun, std::invocable<int64_t> SummatoryFun>
constexpr auto primeSum(size_t limit, Fun f, SummatoryFun F)
{
    return primeSumTable(limit, std::move(f), std::move(F))[limit];
}

/// Calculates `(1 ≤ p ≤ limit, p prime) * p`.
template <typename T = int64_t> constexpr T primeSum(size_t limit)
{
    return primeSum(limit, [](size_t n) -> T { return n; }, [](auto &&n) -> T { return sumId<T>(n); });
}

/// Returns a list of pairs `(exp, c)` indicating that `c` primes have exponent `exp` in the factorization of `n!`.
/// O(n^(3/4)). Sublinear version of `factorFactorial`.
inline std::vector<std::pair<int64_t, int64_t>> factorialExponents(int64_t n)
{
    uint32_t const s = isqrt(n);
    std::vector<std::pair<int64_t, int64_t>> res;
    it::primes(2, s)([&](auto p) { res.emplace_back(factorialValuation(n, p), 1); });
    auto const S = primePiTable(n);
    for (int64_t i = n / (s + 1); i >= 1; --i)
    {
        int64_t const c = S[n / i] - S[n / (i + 1)];
        if (c > 0)
            res.emplace_back(i, c);
    }
    return res;
}
} // namespace euler
