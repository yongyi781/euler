#pragma once

#include "decls.hpp"

inline namespace euler
{
template <typename T> struct euclidean_result_t
{
    T gcd;
    T x;
    T y;
};

/// Runs the extended Euclidean algorithm. Returns the triple `(g, x, y)` such that `g = gcd(a, b) = xm + yn`.
template <typename T> constexpr euclidean_result_t<T> extendedEuclidean(T m, T n)
{
    if (m < 1 || n < 1)
        throw std::domain_error("extended_euclidean: arguments must be strictly positive");

    bool swapped = false;
    if (m < n)
    {
        swapped = true;
        std::swap(m, n);
    }
    T u0 = m;
    T u1 = 1;
    T u2 = 0;
    T v0 = n;
    T v1 = 0;
    T v2 = 1;
    T w0;
    T w1;
    T w2;
    while (v0 > 0)
    {
        T q = u0 / v0;
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

    euclidean_result_t<T> result;
    result.gcd = u0;
    if (!swapped)
    {
        result.x = u1;
        result.y = u2;
    }
    else
    {
        result.x = u2;
        result.y = u1;
    }

    return result;
}

/// @brief Computes modular inverse.
/// @param a A number.
/// @param modulus The modulus.
/// @return number^-1 mod modulus.
template <integral2 Ta, integral2 Tm> constexpr std::common_type_t<Ta, Tm> modInverse(Ta a, Tm modulus)
{
    using Tp = std::common_type_t<Ta, Tm>;
    assert(modulus > 0);
    if (modulus == 1)
        return 0;
    // make sure a < modulus:
    a = mod(a, modulus);
    if (a == 0)
        return 0;
    euclidean_result_t<Tp> u = extendedEuclidean(Tp(a), Tp(modulus));
    if (u.gcd > 1)
        return 0;
    // x might not be in the range 0 < x < m, let's fix that:
    while (u.x <= 0)
        u.x += modulus;
    return u.x;
}

/// Modular exponentiation.  Works for negative exponents.
template <typename T, integral2 Te, typename Tm>
constexpr T powm(const T &base, const Te &exponent, Tm modulus, const T &identity = T(1))
{
    if constexpr (integral2<T>)
    {
        if (exponent < 0)
            return modInverse(pow(T(base % modulus), T(-exponent), identity, mod_multiplies(modulus)), modulus);
        return pow(T(base % modulus), exponent, identity, mod_multiplies(modulus));
    }
    else
    {
        assert(exponent >= 0);
        return pow(T(base % modulus), exponent, identity, mod_multiplies(modulus));
    }
}

/// A version of modular exponentiation that handles overflow correctly.
template <typename T, integral2 Te, integral2 Tm>
constexpr T powmSafe(const T &base, const Te &exponent, Tm modulus, const T &identity = T(1))
{
    assert(modulus > 0);
    if constexpr (integral2<T>)
    {
        if (exponent < 0)
            return modInverse(pow(T(base % modulus), T(-exponent), identity, mod_multiplies_safe(modulus)), modulus);
        return pow(T(base % modulus), exponent, identity, mod_multiplies_safe(modulus));
    }
    else
    {
        assert(exponent >= 0);
        return pow(T(base % modulus), exponent, identity, mod_multiplies_safe(modulus));
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
