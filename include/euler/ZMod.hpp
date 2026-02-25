#pragma once

#include "decls.hpp"
#include "it/base.hpp"
#include "prime.hpp"

namespace euler
{
/// Class for integers modulo a modulus. The modulus must be a compile-time constant.
template <integral2 auto M>
    requires(M > 0 && M <= std::numeric_limits<decltype(M)>::max() / 2)
class ZMod
{
  public:
    using value_type = decltype(M);
    static constexpr value_type modulus = M;
    static constexpr bool is_field = isPrime(M);
    static constexpr bool safe_mul_required = (M > std::numeric_limits<decltype(M)>::max() / M);

    ZMod() = default;

    /// Generic constructor to avoid unintentional narrowing conversion bugs!
    template <integral2 T> constexpr ZMod(T value) : value_(mod(value, M)) {}
    /// Convenience constructor for ZMod<N> whenever M divides N.
    template <integral2 auto N>
        requires(N % M == 0)
    constexpr explicit ZMod(ZMod<N> other) : value_(mod(other.value(), M))
    {
    }

    constexpr explicit operator value_type() const { return value_; }

    std::strong_ordering operator<=>(const ZMod &other) const = default;

    constexpr ZMod &operator+=(const ZMod &other)
    {
        value_ += other.value_;
        if (value_ >= M)
            value_ -= M;
        return *this;
    }

    [[nodiscard]] constexpr friend ZMod operator+(ZMod left, const ZMod &right)
    {
        left += right;
        return left;
    }

    constexpr ZMod &operator-=(const ZMod &other)
    {
        value_ = value_ < other.value_ ? value_ + (M - other.value_) : value_ - other.value_;
        return *this;
    }

    [[nodiscard]] constexpr friend ZMod operator-(ZMod left, const ZMod &right)
    {
        left -= right;
        return left;
    }

    constexpr ZMod &operator*=(const ZMod &other)
    {
        if constexpr (safe_mul_required)
        {
            value_ = mulmod(value_, other.value_, M);
        }
        else
        {
            value_ *= other.value_;
            if (value_ >= M)
                value_ %= M;
        }
        return *this;
    }

    [[nodiscard]] constexpr friend ZMod operator*(ZMod left, const ZMod &right)
    {
        left *= right;
        return left;
    }

    /// Division by -1, 1, 2, 3, 4, 6 are fast. Others divisors use pow or modInverse.
    constexpr ZMod &operator/=(const ZMod &other)
    {
        if (value_ % other.value_ == 0)
        {
            value_ /= other.value_;
            return *this;
        }
        return *this *= ~other;
    }

    /// Division by -1, 1, 2, 3, 4, 6 are fast. Others divisors use pow or modInverse.
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

    constexpr ZMod operator-() const { return {true, value_ == 0 ? 0 : M - value_}; }

    /// Multiplicative inverse. Inverses of -1, 1, 2, 3, 4, 6 are fast, the rest are slow.
    constexpr ZMod operator~() const { return inverse(); }

    template <integral2 T> ZMod &operator<<=(T k)
    {
        if (k < 64)
            *this *= 1_u64 << (int)k;
        else
            *this *= ZMod(2).pow(k);
        return *this;
    }

    template <integral2 T> ZMod operator<<(this ZMod x, T k)
    {
        x <<= k;
        return x;
    }

    template <integral2 T> ZMod &operator>>=(T k)
    {
        if (k == 1)
            *this /= 2;
        else if (k == 2)
            *this /= 4;
        else if (k == 3)
            *this /= 8;
        else
            *this /= ZMod(2).pow(k);
        return *this;
    }

    template <integral2 T> ZMod operator>>(this ZMod x, T k)
    {
        x >>= k;
        return x;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const ZMod &x)
    {
        return o << x.value_;
    }

    /// Returns a mutable reference to the internal value. Handle with care!
    constexpr value_type &value() { return value_; }
    constexpr const value_type &value() const { return value_; }
    constexpr value_type balancedValue() const { return value_ <= M / 2 ? value_ : -(M - value_); }

    /// Returns the multiplicative inverse of this number. Inverses of -1, 1, 2, 3, 4, 6 are fast, the rest are slow.
    constexpr ZMod inverse() const
    {
        if (value_ == 0)
            return 0;
        if (value_ == 1)
            return 1;
        if (value_ == M - 1)
            return M - 1;
        if constexpr (M % 2 != 0)
        {
            if (value_ == 2)
                return {true, (M + 1) / 2};
            if (value_ == 4)
                return {true, M % 4 == 1 ? M - (M - 1) / 4 : (M + 1) / 4};
        }
        if constexpr (M % 3 != 0)
            if (value_ == 3)
                return {true, M % 3 == 1 ? M - (M - 1) / 3 : (M + 1) / 3};
        if constexpr (M % 2 != 0 && M % 3 != 0)
            if (value_ == 6)
                return {true, M % 6 == 1 ? M - (M - 1) / 6 : (M + 1) / 6};
        if constexpr (is_field)
            return pow(M - 2);
        else
            return {true, modInverse(value_, M)};
    }

    /// Returns the modular exponentiation of this number.
    template <integral2 E> constexpr ZMod pow(E exponent) const
    {
        if constexpr (is_field)
            exponent = mod(exponent, M - 1);
        if (exponent == 0 || value_ == 1)
            return 1;
        if (value_ == M - 1)
            return exponent % 2 == 0 ? 1 : -1;
        ZMod x = 1;
        ZMod y = *this;
        if constexpr (std::numeric_limits<E>::is_signed)
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
        if (value_ == 0 || value_ == 1)
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
                return r.value_ <= M / 2 ? r : M - r;
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
                if (i == m.value_)
                    return std::nullopt;
            }
            ZMod b = c.pow(ZMod<M - 1>(2).pow((m - i - 1).value_).value());
            ZMod b2 = b * b;
            r *= b;
            t *= b2;
            c = b2;
            m = i;
        }
        return r.value_ <= M / 2 ? r : M - r;
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
    value_type value_;

    // Shortcut the modulus operation if we know value is between 0 and M-1.
    constexpr ZMod(bool /*unused*/, const decltype(M) &value) : value_(value) {}
};

template <integral2 auto M> constexpr size_t hash_value(const ZMod<M> &n) { return boost::hash_value(n.value()); }
} // namespace euler
