#pragma once

#include "decls.hpp"
#include <boost/multiprecision/gmp.hpp>

inline namespace euler
{
/// Holds a result of the extended Euclidean algorithm, which is a triple (g, s, t).
template <typename T> struct euclidean_result_t
{
    T g, s, t;

    auto operator<=>(const euclidean_result_t &other) const = default;

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const euclidean_result_t &x)
    {
        return o << "{ g: " << x.g << ", s: " << x.s << ", t: " << x.t << " }";
    }
};

/// Runs the extended Euclidean algorithm. Returns the triple `(g, x, y)` such that `g = gcd(a, b) = xm + yn`.
template <integral2 T, integral2 U> constexpr auto xgcd(T m, U n)
{
    using Tp = decltype(auto(boost::multiprecision::detail::evaluate_if_expression(m - n)));
    using std::swap;

    Tp u0 = std::move(m), v0 = std::move(n), u1 = 1, v1 = 0, u2 = 0, v2 = 1, q, w0, w1, w2;
    while (v0 != 0)
    {
        q = u0 / v0;
        w0 = u0 - q * v0;
        w1 = u1 - q * v1;
        w2 = u2 - q * v2;
        swap(u0, v0);
        swap(u1, v1);
        swap(u2, v2);
        swap(v0, w0);
        swap(v1, w1);
        swap(v2, w2);
    }
    if (u0 < 0)
    {
        u0 = -u0;
        u1 = -u1;
        u2 = -u2;
    }
    return euclidean_result_t<Tp>{std::move(u0), std::move(u1), std::move(u2)};
}

/// Extended Euclidean algorithm for `mpz_int`. Returns the triple `(g, x, y)` such that `g = gcd(a, b) = xm + yn`.
inline euclidean_result_t<mpz_int> xgcd(const mpz_int &m, const mpz_int &n)
{
    euclidean_result_t<mpz_int> res;
    mpz_gcdext(res.g.backend().data(), res.s.backend().data(), res.t.backend().data(), m.backend().data(),
               n.backend().data());
    return res;
}

/// Returns the remainder of `a` when divided by `modulus`, in the range [0, modulus).
/// Precondition: `modulus > 0`.
template <integral2 T, integral2 Tm> constexpr auto mod(const T &a, const Tm &modulus)
{
    if (a >= 0)
        return a < modulus ? boost::multiprecision::detail::evaluate_if_expression(a) : a % modulus;
    return boost::multiprecision::detail::evaluate_if_expression(modulus - 1 - (-a - 1) % modulus);
}

/// Non-overflowing modular integer multiplication.
template <integral2 Ta, integral2 Tb, integral2 Tm> constexpr auto modmul(const Ta &a, const Tb &b, const Tm &m)
{
    using T = decltype(auto(boost::multiprecision::detail::evaluate_if_expression(a * b % m)));
    using Td = double_integer_t<T>;
    if constexpr (std::same_as<T, Td>)
    {
        return a * b % m;
    }
    else if constexpr (requires(Td result) { __builtin_mul_overflow(a, b, &result); })
    {
        Td result{};
        __builtin_mul_overflow(a, b, &result);
        return T(result % m);
    }
    else
    {
        return T(Td(a) * Td(b) % m);
    }
}

/// Function object for mod plus.
template <integral2 TMod> struct mod_plus
{
    TMod modulus;

    template <typename T, typename U>
    constexpr auto operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) + std::forward<U>(b)))
    {
        return (std::forward<T>(a) + std::forward<U>(b)) % modulus;
    }
};

/// Function object for mod multiplies.
template <typename TMod> struct mod_multiplies
{
    TMod modulus;

    template <typename T, typename U>
    constexpr auto operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) * std::forward<U>(b)))
    {
        return std::forward<T>(a) * std::forward<U>(b) % modulus;
    }
};

/// This version doesn't overflow.
template <integral2 TMod> struct mod_multiplies_safe
{
    TMod modulus;

    template <typename T, typename U>
    constexpr auto operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) * std::forward<U>(b)))
    {
        return modmul(std::forward<T>(a), std::forward<U>(b), modulus);
    }
};

/// Computes `a^-1 mod m`.
template <integral2 T, integral2 U> constexpr auto modInverse(T a, U m)
{
    using Tp = decltype(auto(boost::multiprecision::detail::evaluate_if_expression(a % m)));
    auto [g, s, _] = xgcd(mod(a, m), m);
    if (g > 1)
        return Tp{};
    // s might not be in the range 0 ≤ s < m, let's fix that:
    if (s < 0 || s >= m)
        s += m;
    return s;
}

/// Specialization of `modInverse` for `mpz_int`.
inline auto modInverse(const mpz_int &a, const mpz_int &modulus)
{
    mpz_int res;
    mpz_invert(res.backend().data(), a.backend().data(), modulus.backend().data());
    return res;
}

namespace detail
{
/// Modular exponentiation. Works for negative exponents.
template <typename Multiply, typename T, integral2 E, typename M>
constexpr auto powm(T base, E exponent, M modulus, T identity = T(1))
{
    assert(modulus > 0);
    if constexpr (integral2<T>)
    {
        // Allow negative exponents in this case
        bool const neg = exponent < 0;
        if (neg)
            exponent = -exponent;
        auto res = pow(T(std::move(base) % modulus), std::move(exponent), std::move(identity), Multiply(modulus));
        if (neg)
            return decltype(res)(modInverse(std::move(res), modulus));
        return res;
    }
    else
    {
        assert(exponent >= 0);
        return pow(T(std::move(base) % modulus), std::move(exponent), std::move(identity), Multiply(modulus));
    }
}
} // namespace detail

/// Modular exponentiation. Works for negative exponents.
template <typename T, integral2 E, typename M> constexpr auto powm(T base, E exponent, M modulus, T identity = T(1))
{
    return detail::powm<mod_multiplies<M>>(std::move(base), std::move(exponent), std::move(modulus),
                                           std::move(identity));
}

/// A version of modular exponentiation that does arithmetic in an integer of twice the size of T, to avoid overflow.
template <typename T, integral2 E, typename M> constexpr auto powmSafe(T base, E exponent, M modulus, T identity = T(1))
{
    return detail::powm<mod_multiplies_safe<M>>(std::move(base), std::move(exponent), std::move(modulus),
                                                std::move(identity));
}

/// Builds the table of inverses modulo a prime.
template <integral2 T>
constexpr std::vector<T> modInverseTable(const T &modulus, size_t limit = std::numeric_limits<size_t>::max())
{
    size_t hardLimit = std::min(limit, (size_t)(modulus - 1));
    std::vector<T> result(hardLimit + 1);
    result[1] = 1;
    for (size_t i = 2; i <= hardLimit; ++i)
        result[i] = result[modulus % i] * (modulus - modulus / i) % modulus;
    return result;
}

/// Takes as input an odd prime p and n < p and returns r such that r * r = n [mod p].
/// @param n A number from 0 to p-1.
/// @param p An odd prime number.
/// @return A square root of n (mod p).
template <integral2 Tn, integral2 Tp> constexpr Tp sqrtModp(Tn n, Tp p)
{
    if (n < 0 || n >= p)
        n = mod(n, p);
    if (n == 0)
        return 0;
    Tp s = 0;
    Tp q = p - 1;
    while ((q & 1) == 0)
    {
        q /= 2;
        ++s;
    }
    if (s == 1)
    {
        // Our modulus is 3 mod 4, there's a fast way to do this!
        Tp r = powm(n, (p + 1) >> 2, p);
        if ((r * r) % p == n)
            return r;
        return Tp(-1);
    }
    // Find the first quadratic non-residue z by brute-force search
    Tp z = 1;
    while (powm(++z, (p - 1) >> 1, p) != p - 1)
    {
    }
    Tp c = powm(z, q, p);
    Tp r = powm(n, (q + 1) >> 1, p);
    Tp t = powm(n, q, p);
    Tp m = s;
    while (t != 1)
    {
        Tp tt = t;
        Tp i = 0;
        while (tt != 1)
        {
            tt = (tt * tt) % p;
            ++i;
            if (i == m)
                return Tp(-1);
        }
        Tp b = powm(c, powm(Tp(2), m - i - 1, p - 1), p);
        Tp b2 = (b * b) % p;
        r = (r * b) % p;
        t = (t * b2) % p;
        c = b2;
        m = i;
    }
    if ((r * r) % p == n)
        return r;
    return Tp(-1);
}

/// Returns the modular inverse of `a` modulo 2^64. Much faster than generic modular inverse algorithm.
/// Precondition: `a` must be odd.
inline u64 modInverse_u64(u64 a)
{
    assert(a % 2 == 1);
    u64 const x0 = (3 * a) ^ 2;
    u64 y = 1 - a * x0;
    u64 const x1 = x0 * (1 + y);
    y *= y;
    u64 const x2 = x1 * (1 + y);
    y *= y;
    u64 const x3 = x2 * (1 + y);
    y *= y;
    u64 const x4 = x3 * (1 + y);
    return x4;
}

/// Returns the `k` least significant bits of `number`.
template <std::integral T> T bit_trunc(T number, int k) { return number & (T(1) << k) - 1; }
} // namespace euler
