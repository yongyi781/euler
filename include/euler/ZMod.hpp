#pragma once

#include "decls.hpp"
#include "it/base.hpp"
#include "prime.hpp"

namespace euler
{
/// Class for integers modulo a modulus. The modulus must be a compile-time constant.
/// @tparam SafeMul whether to use safe multiplication.
template <integral2 auto M, bool SafeMul = (uint128_t)M >= std::numeric_limits<decltype(M)>::max() / M>
    requires(M > 0 && M <= std::numeric_limits<decltype(M)>::max() / 2)
class ZMod
{
    decltype(M) _value;

    // Shortcut the modulus operation if we know value is between 0 and M-1.
    constexpr ZMod(bool /*unused*/, const decltype(M) &value) : _value(value) {}

  public:
    using value_type = decltype(M);
    static constexpr value_type modulus = M;
    static constexpr bool isPrimeModulus = isPrime(M);

    ZMod() = default;

    /// Generic constructor to avoid unintentional narrowing conversion bugs!
    template <integral2 T> constexpr ZMod(T value) : _value(mod(value, M)) {}

    constexpr explicit operator value_type() const { return _value; }

    constexpr ZMod &operator+=(const ZMod &other)
    {
        _value += other._value;
        if (_value >= M)
            _value -= M;
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
        _value = _value < other._value ? _value + (M - other._value) : _value - other._value;
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
            _value = modmul(_value, other._value, M);
        }
        else
        {
            _value *= other._value;
            if (_value >= M)
                _value %= M;
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
        if (_value % other._value == 0)
        {
            _value /= other._value;
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

    constexpr ZMod operator-() const { return {true, _value == 0 ? 0 : M - _value}; }

    /// Multiplicative inverse.
    constexpr ZMod operator~() const { return inverse(); }

    std::strong_ordering operator<=>(const ZMod &other) const = default;

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const ZMod &x)
    {
        return o << x._value;
    }

    constexpr value_type value() const { return _value; }
    constexpr value_type balancedValue() const { return _value <= M / 2 ? _value : -(M - _value); }

    /// Returns the multiplicative inverse of this number.
    constexpr ZMod inverse() const
    {
        if constexpr (isPrimeModulus)
            return pow(M - 2);
        else
            return {true, modInverse(_value, M)};
    }

    /// Returns the modular exponentiation of this number.
    template <integral2 T> constexpr ZMod pow(T exponent) const
    {
        if constexpr (isPrimeModulus)
            exponent = mod(exponent, M - 1);
        if (exponent == 0 || _value == 1)
            return 1;
        if (_value == M - 1)
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
        requires(isPrimeModulus)
    {
        if (_value == 0 || _value == 1)
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
                return r._value <= M / 2 ? r : M - r;
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
                if (i == m._value)
                    return std::nullopt;
            }
            ZMod b = c.pow(ZMod<M - 1>(2).pow((m - i - 1)._value).value());
            ZMod b2 = b * b;
            r *= b;
            t *= b2;
            c = b2;
            m = i;
        }
        return r._value <= M / 2 ? r : M - r;
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
};

template <integral2 auto M> constexpr size_t hash_value(const ZMod<M> &n) { return boost::hash_value(n.value()); }
} // namespace euler

template <euler::integral2 auto M>
class boost::multiprecision::number_category<euler::ZMod<M>>
    : public boost::multiprecision::number_category<decltype(M)>
{
};
