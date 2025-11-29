#pragma once

#include "types.hpp"

namespace euler
{
/// Returns a base raised to an integer power. The type of the base needs a multiplication operation
/// defined on it.
template <typename T, integral2 U, std::invocable<T, T> BinaryOp>
constexpr T pow(T base, U exponent, T identity, BinaryOp op)
{
    if (exponent == 0 || base == identity)
        return identity;
    if (exponent == 1)
        return base;
    if constexpr ((integral2<T> || std::floating_point<T>) && std::is_same_v<BinaryOp, std::multiplies<>> &&
                  boost::multiprecision::is_signed_number<T>::value)
        if (base == -identity)
            return exponent % 2 == 0 ? identity : -identity;
    if constexpr (boost::multiprecision::number_category<T>::value != boost::multiprecision::number_kind_unknown)
    {
        if (exponent < 0)
        {
            base = T(1) / std::move(base);
            exponent = -std::move(exponent);
        }
    }
    else
    {
        assert(exponent >= 0);
    }

    T x = std::move(identity);
    T y = std::move(base);
    while (true)
    {
        if ((exponent & 1) != 0)
        {
            x = op(std::move(x), y);
            if (exponent == 1)
                break;
        }
        exponent >>= 1;
        y = op(y, y);
    }

    return x;
}

/// Returns a base raised to an integer power. The type of the base needs a multiplication operation
/// defined on it. This overrides the standard pow functions.
template <typename T, integral2 U>
constexpr T pow(T base, U exponent)
    requires requires {
        base * base;
        T{1};
    }
{
    return pow(base, exponent, T{1}, std::multiplies{});
}

/// Returns whether `a * b ≤ c`, and always returns false if the multiplication overflows.
template <integral2 T, integral2 U, integral2 V> constexpr bool mulLeq(T a, U b, V c)
{
    if constexpr (std::integral<T> || !std::integral<U>)
    {
        V x{};
#if defined(__clang__) || defined(__GNUC__)
        return !__builtin_mul_overflow(a, b, &x) && x <= c;
#else
        return !std::_Mul_overflow(V(a), V(b), x) && x <= c;
#endif
    }
    else
    {
        return a * b <= c;
    }
}

/// Computes the integral square root of a number.
template <integral2 T> constexpr auto isqrt(const T &n)
{
    if constexpr (!std::integral<T>)
    {
        return sqrt(n);
    }
    else
    {
        // boost::multiprecision::sqrt is constexpr, so take advantage of that in a constant-evaluated context.
        if (std::is_constant_evaluated())
            return boost::multiprecision::sqrt(n);
        // This constant is the first input where floor(sqrt(n)) returns the wrong value.
        if (n < 4'503'599'761'588'224)
            return (T)sqrt(n);
        T x = sqrt(n);
        while (x * x > n)
            --x;
        // double sqrt will underestimate starting from 2^106.
        if constexpr (sizeof(T) > 8)
            while ((x + 1) * (x + 1) <= n)
                ++x;
        return x;
    }
}

/// Computes the integral nth root of a number.
template <integral2 T> constexpr T inth_root(T x, int n)
{
    if (n == 1)
        return x;
    if (n == 2)
        return isqrt(std::move(x));
    if (n == 4)
        return isqrt(isqrt(std::move(x)));
    T s = (T)std::pow((double)x, (1.0 + DBL_EPSILON) / n);
    while (pow(s, n) > x)
        --s;
    return s;
}

/// Finds the largest `e` such that `b^e ≤ n`.
template <integral2 T, integral2 U> constexpr int floor_log(T n, U b)
{
    if (n < b)
        return 0;
    if (n < b * b)
        return 1;
    int e = 1;
    T x = b;
    for (; mulLeq(x, b, n); x *= b, ++e)
        ;
    return e;
}

/// Invoke a callable object, and returns true if the callable returns void.
template <typename Callable, typename... Args>
    requires std::invocable<Callable, Args...>
constexpr bool invokeTrueIfVoid(Callable &&f, Args &&...args) noexcept(std::is_nothrow_invocable_v<Callable, Args...>)
{
    if constexpr (std::is_void_v<std::invoke_result_t<Callable, Args...>>)
    {
        std::forward<Callable>(f)(std::forward<Args>(args)...);
        return true;
    }
    else
    {
        return std::forward<Callable>(f)(std::forward<Args>(args)...);
    }
}

/// Converts any integral type to a string in the specified base.
template <integral2 T> constexpr std::string to_string(T n, int base = 10, bool lowercase = false)
{
    bool const neg = n < 0;
    if (neg)
        n = -n;
    std::string_view const digits =
        lowercase ? R"(0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*()[]{}-+=;:'",./<>?\)"
                  : R"(0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz~!@#$%^&*()[]{}-+=;:'",./<>?\)";
    assert((size_t)base <= digits.size());
    std::string s(1 + std::max(0, floor_log(n, base) + neg), '0');
    if (neg)
        s[0] = '-';
    auto it = s.rbegin();
    while (n)
    {
        *it++ = digits[(int)(n % base)];
        n /= base;
    }
    return s;
}

/// Function object for min.
struct minimum
{
    template <typename T, typename U> constexpr auto operator()(T &&a, U &&b) const
    {
        return std::min(std::forward<T>(a), std::forward<U>(b));
    }
};

/// Function object for max.
struct maximum
{
    template <typename T, typename U> constexpr auto operator()(T &&a, U &&b) const
    {
        return std::max(std::forward<T>(a), std::forward<U>(b));
    }
};

#ifdef _WIN32
/// Sets the output code page used by the console associated with the calling process. A console uses its output code
/// page to translate the character values written by the various output functions into the images displayed in the
/// console window.
extern "C" __declspec(dllimport) int __stdcall SetConsoleOutputCP(unsigned wCodePageID);
#endif

/// Sets the console to UTF-8, so that UTF-8 characters display properly on Windows.
inline void setConsoleToUtf8()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
}
} // namespace euler

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, i128 x)
{
    auto const f = o.flags();
    int base = 10;
    if (f & std::ios_base::hex)
        base = 16;
    else if (f & std::ios_base::oct)
        base = 8;
    return o << euler::to_string(x, base);
}

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, u128 x)
{
    auto const f = o.flags();
    int base = 10;
    if (f & std::ios_base::hex)
        base = 16;
    else if (f & std::ios_base::oct)
        base = 8;
    return o << euler::to_string(x, base);
}
