#pragma once

#include "PF.hpp"
#include "it/primes.hpp"
#include "modular_arithmetic.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

inline namespace euler
{
/// Gaussian integer type.
template <integral2 T> struct GI
{
    T re, im;

    GI(T real = {}, T imag = {}) : re(std::move(real)), im(std::move(imag)) {}

    template <integral2 U> GI<T> operator+(const GI<U> &other) const { return {re + T(other.re), im + T(other.im)}; }
    template <integral2 U> GI<T> &operator+=(const GI<U> &other) { return *this = *this + other; }
    template <integral2 U> GI<T> operator-(const GI<U> &other) const { return {re - T(other.re), im - T(other.im)}; }
    template <integral2 U> GI<T> &operator-=(const GI<U> &other) { return *this = *this - other; }
    template <integral2 U> GI<T> operator*(const GI<U> &other) const
    {
        return {re * T(other.re) - im * T(other.im), re * T(other.im) + im * T(other.re)};
    }
    template <integral2 U> GI<T> &operator*=(const GI<U> &other) { return *this = *this * other; }
    /// Scalar multiplication
    template <integral2 U> GI<T> operator*(const U &scalar) const { return {re * scalar, im * scalar}; }
    template <integral2 U> GI<T> &operator*=(const U &scalar) { return *this = *this * scalar; }
    /// Scalar division (floor division)
    template <integral2 U> GI<T> operator/(const U &scalar) const { return {re / scalar, im / scalar}; }
    template <integral2 U> GI<T> &operator/=(const U &scalar) { return *this = *this / scalar; }
    GI operator-() const { return {-re, -im}; }

    // Equality operator
    auto operator<=>(const GI &other) const = default;

    // Norm and absolute value
    [[nodiscard]] T norm() const { return re * re + im * im; }
    [[nodiscard]] double abs() const { return std::sqrt((double)norm()); }

    // Conjugate
    [[nodiscard]] GI conj() const { return {re, -im}; }

    // Division (returns a pair of GaussianIntegers for quotient and remainder)
    [[nodiscard]] std::pair<GI, GI> div(const GI &other) const
    {
        if (other == GI{})
            throw std::overflow_error("Gaussian integer division by zero");

        // (a+bi)/(c+di) = (ac+bd)/(c^2+d^2) + i(bc-ad)/(c^2+d^2)
        T denominator = other.norm();
        T num_re = re * other.re + im * other.im;
        T num_im = im * other.re - re * other.im;

        T q_re = num_re / denominator;
        T q_im = num_im / denominator;

        GI quotient = {q_re, q_im};
        GI remainder = *this - quotient * other;

        return {std::move(quotient), std::move(remainder)};
    }

    // Stream insertion operator for printing
    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const GI &z)
    {
        if (z.re != 0)
            o << z.re;
        if (z.im != 0)
        {
            if (z.re && z.im > 0)
                o << '+';
            o << z.im << "i";
        }
        return o;
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

/// Enumerates Gaussian integers with a given norm.
template <integral2 T = i64> class GINorms
{
    std::vector<GI<T>> prime_solutions;

  public:
    using value_type = GI<T>;

    GINorms(size_t prime_limit) : prime_solutions(prime_limit)
    {
        auto const ps = it::primes(prime_limit).filter([](i64 p) { return p % 4 != 3; }).to();
        std::for_each(std::execution::par, ps.begin(), ps.end(),
                      [&](i64 p) { prime_solutions[p] = primeNormGI((T)p); });
    }

    /// Enumerates Gaussian integers with a given norm, using the given prime factorization.
    template <integral2 Z, typename Fun> void enumerate(const PF<Z> &fac, Fun callback) const
    {
        // Check early termination condition
        if (std::ranges::any_of(fac, [](auto &&pe) { return pe.first % 4 == 3 && pe.second % 2 != 0; }))
            return;
        auto const loop = [&](this auto &&self, auto it, GI<T> z) -> void {
            if (it == fac.end())
            {
                z = rep(z);
                if (z.re >= z.im)
                    callback(z);
                return;
            }
            auto const [p, e] = *it;
            if (p == 2)
                self(it + 1, z * pow(GI<T>{1, 1}, e));
            else if (p % 4 == 3)
            {
                if (e % 2 == 0)
                    self(it + 1, z * pow(p, e / 2));
            }
            else
            {
                auto const w = prime_solutions[p];
                for (int i = 0; i <= e; ++i)
                    self(it + 1, z * pow(w, e - i) * pow(w.conj(), i));
            }
        };
        loop(fac.begin(), {1, 0});
    }

  private:
    /// Returns a multiple-of-90 degree rotation of z with positive coordinates.
    static GI<T> rep(GI<T> z)
    {
        if (z.re > 0 && z.im >= 0)
            return z;
        if (z.im > 0 && z.re <= 0)
            return {z.im, -z.re};
        if (z.re < 0 && z.im <= 0)
            return -z;
        return {-z.im, z.re};
    }
};
} // namespace euler
