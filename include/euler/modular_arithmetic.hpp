#pragma once

#include "decls.hpp"
#include <boost/multiprecision/gmp.hpp>

namespace euler
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

    Tp g0 = std::move(m), g1 = std::move(n), s0 = 1, s1 = 0, t0 = 0, t1 = 1, q, g2, s2, t2;
    while (g1 != 0)
    {
        q = g0 / g1;
        g2 = g0 - q * g1;
        s2 = s0 - q * s1;
        t2 = t0 - q * t1;
        swap(g0, g1);
        swap(s0, s1);
        swap(t0, t1);
        swap(g1, g2);
        swap(s1, s2);
        swap(t1, t2);
    }
    if constexpr (boost::multiprecision::is_signed_number<Tp>::value)
        if (g0 < 0)
        {
            g0 = -g0;
            s0 = -s0;
            t0 = -t0;
        }
    return euclidean_result_t<Tp>{std::move(g0), std::move(s0), std::move(t0)};
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
    if constexpr (boost::multiprecision::is_unsigned_number<T>::value)
    {
        return a < modulus ? boost::multiprecision::detail::evaluate_if_expression(a) : a % modulus;
    }
    else
    {
        if (a >= 0)
            return a < modulus ? boost::multiprecision::detail::evaluate_if_expression(a) : a % modulus;
        return boost::multiprecision::detail::evaluate_if_expression(modulus - 1 - (-a - 1) % modulus);
    }
}

/// Non-overflowing modular integer multiplication.
template <integral2 Ta, integral2 Tb, integral2 Tm> constexpr auto mulmod(const Ta &a, const Tb &b, const Tm &m)
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
        return mulmod(std::forward<T>(a), std::forward<U>(b), modulus);
    }
};

/// Computes `a^-1 mod m`.
template <integral2 T, integral2 U> constexpr auto modInverse(T a, U m)
{
    using Tp = decltype(auto(boost::multiprecision::detail::evaluate_if_expression(a % m)));
    auto const [g, s, _] = xgcd(mod(a, m), m);
    if (g != 1)
        return Tp{};
    return s + m < m ? s + m : s;
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
    if constexpr (integral2<T>)
    {
        bool neg = false;
        if constexpr (boost::multiprecision::is_signed_number<E>::value)
        {
            neg = exponent < 0;
            if (neg)
                exponent = -exponent;
        }
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

/// Builds the table of inverses in a type T (typically `ZMod<M>`).
template <typename T> std::vector<T> inverseTable(size_t limit)
{
    std::vector<T> res(limit + 1);
    res[1] = 1;
    tbb::parallel_for(
        tbb::blocked_range(2UZ, limit + 1, (1UZ << 17) / sizeof(T)),
        [&](auto r) {
            size_t const start = r.begin();
            thread_local std::vector<T> prefix;
            prefix.resize(r.size() + 1);
            prefix[0] = 1;
            T prod = 1;
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                prod *= i;
                prefix[i - start + 1] = prod;
            }
            T inv = T(1) / prod;
            for (size_t i = r.end(); i-- != r.begin();)
            {
                res[i] = inv * prefix[i - start];
                inv *= i;
            }
        },
        tbb::simple_partitioner{});
    return res;
}

/// Builds the table of inverses modulo a prime.
template <integral2 T> constexpr std::vector<T> modInverseTable(const T &modulus, size_t limit)
{
    std::vector<T> result(limit + 1);
    result[1] = 1;
    for (size_t i = 2; i <= limit; ++i)
        result[i] = result[modulus % i] * (modulus - modulus / i) % modulus;
    return result;
}

/// Returns the smaller solution to `√n` (mod p), or 0 if there is no square root. Uses the Tonelli-Shanks algorithm.
/// Requires `p` to be prime.
template <integral2 T, integral2 U> constexpr auto sqrtModp(T n, U p)
{
    using Tp = std::common_type_t<T, U>;
    Tp x = n;
    if (x < 0 || n >= p)
        x = mod(x, p);
    if (x <= 1)
        return x;
    int s = 0;
    Tp q = p - 1;
    for (; (q & 1) == 0; ++s, q >>= 1)
        ;
    if (s == 1)
    {
        Tp res = euler::powm(x, (p + 1) >> 2, p);
        return res * res % p == x ? std::min(res, p - res) : Tp(0);
    }
    if (s == 2)
    {
        Tp res = euler::powm(x, (p + 3) >> 3, p);
        if (res * res % p == x)
            return std::min(res, p - res);
        res = res * euler::powm(Tp(2), (p - 1) >> 2, p) % p;
        return res * res % p == x ? std::min(res, p - res) : Tp(0);
    }
    if (euler::powm(x, (p - 1) >> 1, p) != 1) // Perform Legendre check now
        return Tp(0);
    Tp z = 2;
    for (; euler::powm(z, (p - 1) >> 1, p) != p - 1; ++z)
        ;
    Tp res = euler::powm(x, (q + 1) >> 1, p);
    Tp c = euler::powm(z, q, p);
    Tp t = euler::powm(x, q, p);
    int m = s;
    while (t != 1)
    {
        Tp tt = t;
        int i = 0;
        for (; tt != 1; ++i, tt = tt * tt % p)
            ;
        Tp b = c;
        for (int k = 0; k < m - i - 1; ++k)
            b = b * b % p;
        res = res * b % p;
        c = b * b % p;
        t = t * c % p;
        m = i;
    }
    return std::min(res, p - res);
}

/// Computes the multiplicative modular inverse of `n` modulo `2^k`, where `k` is the bit width of T (e.g., `2^64` for
/// `u64`). Precondition: `n` must be odd.
template <std::unsigned_integral T> constexpr T bitInverse(T n)
{
    assert(n % 2 == 1);
    T x = (3 * n) ^ 2;
    T y = 1 - n * x;
    x *= (1 + y);
    y *= y;
    x *= (1 + y);
    y *= y;
    x *= (1 + y);
    if constexpr (sizeof(T) > 4)
    {
        y *= y;
        x *= (1 + y);
    }
    if constexpr (sizeof(T) > 8)
    {
        y *= y;
        x *= (1 + y);
    }
    return x;
}

/// Returns the `k` least significant bits of `number`.
template <std::integral T> constexpr T bit_trunc(T number, int k) { return number & (T(1) << k) - 1; }
} // namespace euler
