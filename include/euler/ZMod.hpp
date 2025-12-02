#pragma once

#include "decls.hpp"
#include "it/base.hpp"
#include "prime.hpp"

namespace euler
{
/// Class for integers modulo a modulus. The modulus must be a compile-time constant.
/// @tparam SafeMul whether to use safe multiplication.
template <integral2 auto M, bool SafeMul = (M > std::numeric_limits<decltype(M)>::max() / M)>
    requires(M > 0 && M <= std::numeric_limits<decltype(M)>::max() / 2)
class ZMod
{
  public:
    using value_type = decltype(M);
    static constexpr value_type modulus = M;
    static constexpr bool is_field = isPrime(M);

    ZMod() = default;

    /// Generic constructor to avoid unintentional narrowing conversion bugs!
    template <integral2 T> constexpr ZMod(T value) : m_value(mod(value, M)) {}
    /// Convenience constructor for ZMod<N> whenever M divides N.
    template <integral2 auto N>
        requires(N % M == 0)
    constexpr ZMod(ZMod<N> other) : m_value(mod(other.value(), M))
    {
    }

    constexpr explicit operator value_type() const { return m_value; }

    constexpr ZMod &operator+=(const ZMod &other)
    {
        m_value += other.m_value;
        if (m_value >= M)
            m_value -= M;
        return *this;
    }

    [[nodiscard]] constexpr friend ZMod operator+(ZMod left, const ZMod &right)
    {
        left += right;
        return left;
    }

    // constexpr friend ZMod operator+(const integral2 auto &left, const ZMod &right) { return ZMod{left} + right; }

    constexpr ZMod &operator-=(const ZMod &other)
    {
        m_value = m_value < other.m_value ? m_value + (M - other.m_value) : m_value - other.m_value;
        return *this;
    }

    [[nodiscard]] constexpr friend ZMod operator-(ZMod left, const ZMod &right)
    {
        left -= right;
        return left;
    }

    constexpr ZMod &operator*=(const ZMod &other)
    {
        if constexpr (SafeMul)
        {
            m_value = mulmod(m_value, other.m_value, M);
        }
        else
        {
            m_value *= other.m_value;
            if (m_value >= M)
                m_value %= M;
        }
        return *this;
    }

    [[nodiscard]] constexpr friend ZMod operator*(ZMod left, const ZMod &right)
    {
        left *= right;
        return left;
    }

    constexpr ZMod &operator/=(const ZMod &other)
    {
        if (m_value % other.m_value == 0)
        {
            m_value /= other.m_value;
            return *this;
        }
        return *this *= ~other;
    }

    [[nodiscard]] constexpr friend ZMod operator/(ZMod left, const ZMod &right)
    {
        left /= right;
        return left;
    }

    constexpr ZMod &operator++() { return *this += 1; }

    constexpr ZMod operator++(int)
    {
        ZMod result = *this;
        ++*this;
        return result;
    }

    constexpr ZMod &operator--() { return *this -= 1; }

    constexpr ZMod operator--(int)
    {
        ZMod result = *this;
        ++*this;
        return result;
    }

    constexpr ZMod operator-() const { return {true, m_value == 0 ? 0 : M - m_value}; }

    /// Multiplicative inverse.
    constexpr ZMod operator~() const { return inverse(); }

    std::strong_ordering operator<=>(const ZMod &other) const = default;

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const ZMod &x)
    {
        return o << x.m_value;
    }

    /// Returns a mutable reference to the internal value. Handle with care!
    constexpr value_type &value() { return m_value; }
    constexpr const value_type &value() const { return m_value; }
    constexpr value_type balancedValue() const { return m_value <= M / 2 ? m_value : -(M - m_value); }

    /// Returns the multiplicative inverse of this number.
    constexpr ZMod inverse() const
    {
        if constexpr (is_field)
            return pow(M - 2);
        else
            return {true, modInverse(m_value, M)};
    }

    /// Returns the modular exponentiation of this number.
    template <integral2 T> constexpr ZMod pow(T exponent) const
    {
        if constexpr (is_field)
            exponent = mod(exponent, M - 1);
        if (exponent == 0 || m_value == 1)
            return 1;
        if (m_value == M - 1)
            return exponent % 2 == 0 ? 1 : -1;
        ZMod x = 1;
        ZMod y = *this;
        if (exponent < 0)
        {
            y = ~y;
            exponent = -exponent;
        }
        while (true)
        {
            if (exponent & 1)
            {
                x *= y;
                if (exponent == 1)
                    break;
            }
            exponent >>= 1;
            y *= y;
        }
        return x;
    }

    /// Returns the smaller square root of this residue class if there is one, or nullopt if there
    /// is none. Assumes M is prime.
    constexpr std::optional<ZMod> sqrt() const
        requires(is_field)
    {
        if (m_value == 0 || m_value == 1)
            return *this;
        value_type s = 0;
        value_type q = M - 1;
        while ((q & 1) == 0)
        {
            q /= 2;
            ++s;
        }
        if (s == 1)
        {
            // Our modulus is 3 mod 4, there's a fast way to do this!
            ZMod r = pow((M + 1) / 4);
            if (r * r == *this)
                return r.m_value <= M / 2 ? r : M - r;
            return std::nullopt;
        }
        // Find the first quadratic non-residue z by brute-force search
        ZMod z = 1;
        while ((++z).pow((M - 1) / 2) != M - 1)
        {
        }
        ZMod c = z.pow(q);
        ZMod r = pow((q + 1) / 2);
        ZMod t = pow(q);
        ZMod m = s;
        while (t != 1)
        {
            ZMod tt = t;
            value_type i = 0;
            while (tt != 1)
            {
                tt *= tt;
                ++i;
                if (i == m.m_value)
                    return std::nullopt;
            }
            ZMod b = c.pow(ZMod<M - 1>(2).pow((m - i - 1).m_value).value());
            ZMod b2 = b * b;
            r *= b;
            t *= b2;
            c = b2;
            m = i;
        }
        return r.m_value <= M / 2 ? r : M - r;
    }

    /// Calculates factorial of n mod M. Does not use Wilson's theorem.
    template <execution_policy Exec, integral2 T> static ZMod factorial(Exec &&exec, T n)
    {
        return it::range(T(1), n).map([](T k) { return ZMod{k}; }).product(std::forward<Exec>(exec));
    }

    /// Calculates factorial of n mod M. Does not use Wilson's theorem.
    template <integral2 T> static constexpr ZMod factorial(T n)
    {
        return it::range(T(1), n).map([](T k) { return ZMod{k}; }).product();
    }

    /// Calculates binomial coefficient (n choose r) mod M. Does not use Lucas's theorem.
    template <execution_policy Exec, integral2 T> static ZMod binomial(Exec &&exec, T n, T r)
    {
        if (r <= 0 || n == 0)
            return ZMod(r == 0);
        if (n < 0)
            return euler::pow(-1, r) * binomial(exec, r - n - 1, r);
        if (r > n)
            return 0;
        r = std::min(r, n - r);
        return it::range(T(0), r - 1).map([&](auto k) { return ZMod{n - k}; }).product(std::forward<Exec>(exec)) /
               it::range(T(0), r - 1).map([&](auto k) { return ZMod{k + 1}; }).product(std::forward<Exec>(exec));
    }

    /// Calculates binomial coefficient (n choose r) mod M. Does not use Lucas's theorem.
    template <integral2 T> static constexpr ZMod binomial(T n, T r)
    {
        if (r <= 0 || n == 0)
            return ZMod(r == 0);
        if (n < 0)
            return euler::pow(-1, r) * binomial(r - n - 1, r);
        if (r > n)
            return 0;
        r = std::min(r, n - r);
        return it::range(T(0), r - 1).map([&](auto k) { return ZMod{n - k}; }).product() /
               it::range(T(0), r - 1).map([&](auto k) { return ZMod{k + 1}; }).product();
    }

  private:
    value_type m_value;

    // Shortcut the modulus operation if we know value is between 0 and M-1.
    constexpr ZMod(bool /*unused*/, const decltype(M) &value) : m_value(value) {}
};

template <integral2 auto M> constexpr size_t hash_value(const ZMod<M> &n) { return boost::hash_value(n.value()); }
} // namespace euler
