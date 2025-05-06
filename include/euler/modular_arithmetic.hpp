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
/// Precondition: T is signed.
template <typename T> constexpr euclidean_result_t<T> xgcd(T m, T n)
{
    T u0 = std::move(m), v0 = std::move(n), u1 = 1, v1 = 0, u2 = 0, v2 = 1, q, w0, w1, w2;
    while (v0 != 0)
    {
        q = u0 / v0;
        w0 = u0 - q * v0;
        w1 = u1 - q * v1;
        w2 = u2 - q * v2;
        u0 = v0;
        u1 = v1;
        u2 = v2;
        v0 = w0;
        v1 = w1;
        v2 = w2;
    }
    if (u0 < 0)
    {
        u1 = -u1;
        u2 = -u2;
    }
    return {std::move(u0), std::move(u1), std::move(u2)};
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

/// Computes modular inverse.
/// @param a A number.
/// @param modulus The modulus.
/// @return a^-1 mod modulus.
template <integral2 Ta, integral2 Tm> constexpr auto modInverse(Ta a, Tm modulus)
{
    using T = decltype(auto(boost::multiprecision::detail::evaluate_if_expression(a % modulus)));
    using Ts = std::make_signed_t<T>;
    assert(modulus > 0);
    if (modulus == 1)
        return Ts{};
    // make sure a < modulus:
    Ts a_small = Ts(mod(a, modulus));
    if (a_small == 0)
        return Ts{};
    auto [g, s, _] = xgcd(a_small, Ts(modulus));
    if (g > 1)
        return Ts{};
    // s might not be in the range 0 ≤ s < m, let's fix that:
    if (s < 0)
        s += modulus;
    return s;
}

/// Specialization of `modInverse` for `mpz_int`.
inline auto modInverse(const mpz_int &a, const mpz_int &modulus)
{
    mpz_int res;
    mpz_invert(res.backend().data(), a.backend().data(), modulus.backend().data());
    return res;
}

/// Modular exponentiation. Works for negative exponents.
template <typename T, integral2 Te, typename Tm>
constexpr auto powm(T base, Te exponent, const Tm &modulus, T identity = T(1))
{
    assert(modulus > 0);
    if constexpr (integral2<T>)
    {
        if (exponent < 0)
            return modInverse(pow(boost::multiprecision::detail::evaluate_if_expression(std::move(base) % modulus),
                                  boost::multiprecision::detail::evaluate_if_expression(-std::move(exponent)),
                                  std::move(identity), mod_multiplies(modulus)),
                              modulus);
        return pow(boost::multiprecision::detail::evaluate_if_expression(std::move(base) % modulus),
                   std::move(exponent), std::move(identity), mod_multiplies(modulus));
    }
    else
    {
        assert(exponent >= 0);
        return pow(boost::multiprecision::detail::evaluate_if_expression(std::move(base) % modulus),
                   std::move(exponent), std::move(identity), mod_multiplies(modulus));
    }
}

/// A version of modular exponentiation that handles overflow correctly.
template <typename T, integral2 Te, integral2 Tm>
constexpr auto powmSafe(T base, Te exponent, const Tm &modulus, const T &identity = T(1))
{
    assert(modulus > 0);
    if constexpr (integral2<T>)
    {
        if (exponent < 0)
            return modInverse(pow(boost::multiprecision::detail::evaluate_if_expression(std::move(base) % modulus),
                                  boost::multiprecision::detail::evaluate_if_expression(-std::move(exponent)),
                                  std::move(identity), mod_multiplies_safe(modulus)),
                              modulus);
        return pow(boost::multiprecision::detail::evaluate_if_expression(std::move(base) % modulus),
                   std::move(exponent), std::move(identity), mod_multiplies_safe(modulus));
    }
    else
    {
        assert(exponent >= 0);
        return pow(boost::multiprecision::detail::evaluate_if_expression(std::move(base) % modulus),
                   std::move(exponent), std::move(identity), mod_multiplies_safe(modulus));
    }
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
        Tp r = powm(n, (p + 1) / 4, p);
        if ((r * r) % p == n)
            return r;
        return Tp(-1);
    }
    // Find the first quadratic non-residue z by brute-force search
    Tp z = 1;
    while (powm(++z, (p - 1) / 2, p) != p - 1)
    {
    }
    Tp c = powm(z, q, p);
    Tp r = powm(n, (q + 1) / 2, p);
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
        Tp b = powm(c, powm(2LL, m - i - 1, p - 1), p);
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
} // namespace euler
