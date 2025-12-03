#pragma once

#include "it/primes.hpp"
#include "math.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace euler
{
/// Gaussian integer type.
template <integral2 T> struct GI
{
    T re, im;

    constexpr GI(T real = {}, T imag = {}) : re(std::move(real)), im(std::move(imag)) {}

    template <integral2 U> constexpr friend auto operator+(const GI<T> &left, const GI<U> &right)
    {
        using V = std::common_type_t<T, U>;
        return GI<V>{V(left.re + right.re), V(left.im + right.im)};
    }
    template <integral2 U> constexpr GI<T> &operator+=(const GI<U> &other) { return *this = *this + other; }

    template <integral2 U> constexpr friend auto operator-(const GI<T> &left, const GI<U> &right)
    {
        using V = std::common_type_t<T, U>;
        return GI<V>{V(left.re - right.re), V(left.im - right.im)};
    }
    template <integral2 U> constexpr GI<T> &operator-=(const GI<U> &other) { return *this = *this - other; }
    template <integral2 U> constexpr GI<T> operator*(const GI<U> &other) const
    {
        return {re * T(other.re) - im * T(other.im), re * T(other.im) + im * T(other.re)};
    }
    template <integral2 U> constexpr GI<T> &operator*=(const GI<U> &other) { return *this = *this * other; }

    /// Truncated division in the Gaussian integers.
    template <integral2 U> constexpr GI<T> operator/(const GI<U> &other) const { return div(other).first; }
    template <integral2 U> constexpr GI<T> &operator/=(const GI<U> &other) { return *this = *this / other; }
    template <integral2 U> constexpr GI<T> operator%(const GI<U> &other) const { return div(other).second; }
    template <integral2 U> constexpr GI<T> &operator%=(const GI<U> &other) { return *this = *this % other; }

    /// Scalar multiplication
    template <integral2 U> constexpr GI<T> operator*(const U &scalar) const { return {re * scalar, im * scalar}; }
    template <integral2 U> constexpr GI<T> &operator*=(const U &scalar) { return *this = *this * scalar; }
    template <integral2 U> friend constexpr GI operator*(const U &scalar, const GI &z) { return z * scalar; }
    /// Scalar division (floor division)
    template <integral2 U> constexpr GI<T> operator/(const U &scalar) const { return {re / scalar, im / scalar}; }
    template <integral2 U> constexpr GI<T> &operator/=(const U &scalar) { return *this = *this / scalar; }
    constexpr GI operator-() const { return {-re, -im}; }

    // Equality operator
    constexpr auto operator<=>(const GI &other) const = default;

    // Norm and absolute value
    [[nodiscard]] constexpr T norm() const { return re * re + im * im; }
    [[nodiscard]] double abs() const { return std::sqrt((double)norm()); }

    // Conjugate
    [[nodiscard]] constexpr GI conj() const { return {re, -im}; }

    /// Rotate by a power of i so that (re > 0 and im ≥ 0).
    [[nodiscard]] constexpr GI canonicalAssociate() const
    {
        if (re > 0 && im >= 0)
            return *this;
        if (im > 0 && re <= 0)
            return {im, -re};
        if (re < 0 && im <= 0)
            return -*this;
        return {-im, re};
    }

    // Division with remainder. This is normalized so that the remainder has the smallest norm among possible
    // remainders.
    template <integral2 U> [[nodiscard]] constexpr std::pair<GI, GI> div(const GI<U> &other) const
    {
        T const d = other.norm();
        if (d == 0)
            throw std::overflow_error("Gaussian integer division by zero");
        auto const [a, b] = *this * other.conj();

        GI q{floorDiv(a + (d - 1) / 2, d), floorDiv(b + (d - 1) / 2, d)};
        GI r = *this - q * other;

        return {std::move(q), std::move(r)};
    }

    // Stream insertion operator for printing
    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const GI &z)
    {
        if (z.re != 0 || z.im == 0)
            o << z.re;
        if (z.im != 0)
        {
            if (z.re != 0 && z.im > 0)
                o << '+';
            if (z.im == -1)
                o << "-i";
            else if (z.im == 1)
                o << "i";
            else
                o << z.im << "i";
        }
        return o;
    }

    friend size_t hash_value(const GI &z)
    {
        size_t seed = 0;
        boost::hash_combine(seed, z.re);
        boost::hash_combine(seed, z.im);
        return seed;
    }
};

/// Returns the unique Gaussian integer a + bi with norm p such that a ≥ b. Precondition: p must be prime.
template <integral2 T> inline GI<T> primeNormGI(T p)
{
    if (p < 5'000'000)
    {
        // Brute force
        for (T a = 0; a * a < p; ++a)
        {
            T x = isqrt(p - a * a);
            if (x * x + a * a == p)
                return {x, a};
        }
        throw std::runtime_error("No solution found");
    }
    // Cornacchia’s algorithm for p > 5'000'000
    T const t = sqrtModp(p - 1, p);
    T x = p, a = t;
    while (a * a > p)
        std::tie(a, x) = std::pair<T, T>{x % a, a};
    T b = isqrt(p - a * a);
    if (a < b)
        std::swap(a, b);
    return {a, b};
}

/// Enumerates Gaussian integers with a given norm, using the given prime factorization.
/// If `includeAll` is true, then both `(x, y)` and `(y, x)` will be reported. Otherwise, only `(x, y)` with `x ≥ y`
/// will be reported.
template <integral2 T = i64, integral2 Z, typename Fun>
void enumGIWithNorm(const PF<Z> &fac, Fun callback, bool includeAll = false)
{
    // Check early termination condition
    if (std::ranges::any_of(fac, [](auto &&pe) { return pe.first % 4 == 3 && pe.second % 2 != 0; }))
        return;

    auto const v = mapv(fac, [&](auto &&pe) {
        auto &&[p, e] = pe;
        return std::tuple{p, e, p % 4 == 3 ? GI<T>{} : primeNormGI((T)p)};
    });
    auto const loop = [&](this auto &&self, auto it, GI<T> z) -> void {
        if (it == v.end())
        {
            auto const a = z.canonicalAssociate();
            if (includeAll || a.re >= a.im)
                callback(a);
            return;
        }
        auto const [p, e, w] = *it;
        if (p == 2)
            self(it + 1, z * pow(w, e));
        else if (p % 4 == 3)
            self(it + 1, z * pow(p, e / 2));
        else
        {
            GI<T> q = 1;
            for (int i = 0; i <= e; ++i, q *= w.conj())
                self(it + 1, z * pow(w, e - i) * q);
        }
    };
    loop(v.begin(), {1, 0});
}

/// Enumerates Gaussian integers with a given norm. Use this if you want to precompute prime solutions up front.
template <integral2 T = i64> class GINorms
{
    std::vector<GI<T>> prime_solutions;

  public:
    using value_type = GI<T>;

    GINorms(size_t prime_limit) : prime_solutions(prime_limit)
    {
        auto const ps = it::primes(2, prime_limit).filter([](i64 p) { return p % 4 != 3; }).to();
        std::for_each(std::execution::par, ps.begin(), ps.end(),
                      [&](i64 p) { prime_solutions[p] = primeNormGI((T)p); });
    }

    /// Enumerates Gaussian integers with a given norm, using the given prime factorization.
    template <integral2 Z, typename Fun> void enumerate(const PF<Z> &fac, Fun callback, bool includeAll = false) const
    {
        // Check early termination condition
        if (std::ranges::any_of(fac, [](auto &&pe) { return pe.first % 4 == 3 && pe.second % 2 != 0; }))
            return;
        auto const loop = [&, the = this](this auto &&self, auto it, GI<T> z) -> void {
            if (it == fac.end())
            {
                auto const a = z.canonicalAssociate();
                if (includeAll || a.re >= a.im)
                    callback(a);
                return;
            }
            auto const [p, e] = *it;
            auto const w = the->prime_solutions[p];
            if (p == 2)
                self(it + 1, z * pow(w, e));
            else if (p % 4 == 3)
                self(it + 1, z * pow(p, e / 2));
            else
            {
                GI<T> q = 1;
                for (int i = 0; i <= e; ++i, q *= w.conj())
                    self(it + 1, z * pow(w, e - i) * q);
            }
        };
        loop(fac.begin(), {1, 0});
    }
};
} // namespace euler
