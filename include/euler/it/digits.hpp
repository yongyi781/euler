#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates the digits of a number in a specified base, which must be ≥ 2.
template <integral2 T = int64_t, integral2 TBase = int> class digits : public it_base
{
    T _n;
    TBase _base;

  public:
    using value_type = TBase;

    digits() = default;
    constexpr digits(T n, TBase base = 10) : _n(n), _base(base) { assert(base >= 2); }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if constexpr (boost::integer_traits<T>::digits > 128)
            if (_base == 10)
            {
                for (char c : to_string(_n) | std::views::reverse)
                    if (!callbackResult(f, TBase(c - '0')))
                        return result_break;
                return result_continue;
            }
        T n = _n;
        while (n > 0)
        {
            if (!callbackResult(f, TBase(n % _base)))
                return result_break;
            n /= _base;
        }
        return result_continue;
    }
};
} // namespace it

template <std::endian Endian = std::endian::little, integral2 TBase = int>
constexpr std::vector<TBase> digits(const integral2 auto &num, TBase base = 10)
{
    std::vector<TBase> result;
    result.reserve(8); // Easy optimization
    it::digits(num, base)([&](TBase d) { result.push_back(d); });
    if constexpr (Endian == std::endian::big)
        std::ranges::reverse(result);
    return result;
}

template <integral2 T = i64, std::endian Endian = std::endian::little, std::ranges::range Range>
constexpr T fromDigits(Range &&digits, int base = 10)
{
    T res = 0;
    if constexpr (Endian == std::endian::big)
    {
        for (auto &&d : digits)
            res = res * base + d;
    }
    else
    {
        T p = 1;
        for (auto &&d : digits)
        {
            res += p * d;
            p *= base;
        }
    }
    return res;
}

/// Shortcut for `it::digits(num, base).size()`.
template <integral2 T, integral2 TBase = int> constexpr size_t countDigits(const T &num, TBase base = 10)
{
    return it::digits(num, base).size();
}

/// Sums the digits of a number.
template <integral2 T, integral2 TBase = int> constexpr TBase sumDigits(const T &num, TBase base = 10)
{
    return it::digits(num, base).sum();
}

/// Returns the number with the reversed digits of n.
template <integral2 T, integral2 TBase = int> constexpr T reverseDigits(const T &n, TBase base = 10)
{
    return it::digits(n, base).reduce(T(0), [&](auto &&a, auto &&b) -> T { return base * a + b; });
}

/// Returns whether a number is a palindrome.
template <integral2 T, integral2 TBase = int> constexpr bool isPalindrome(const T &n, TBase base = 10)
{
    return n == reverseDigits(n, base);
}

/// Calculates the digit signature of a given number.
template <int Base = 10> constexpr std::array<int, Base> getDigitSignature(const integral2 auto &n)
{
    std::array<int, Base> result{};
    it::digits(n, Base)([&](int d) { ++result[d]; });
    return result;
}
} // namespace euler
