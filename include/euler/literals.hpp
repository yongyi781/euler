#pragma once

#include "types.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <cstdint>
#include <stdexcept>

namespace euler
{
/// Converts a string to a number in a constant-evaluated context.
template <typename T> consteval T fromString(const char *str)
{
    T res{};
    int base = 10;

    // Detect base from prefix if present.
    if (str[0] == '0')
    {
        if (str[1] == 'x' || str[1] == 'X')
        {
            base = 16;
            str += 2;
        }
        else if (str[1] == 'b' || str[1] == 'B')
        {
            base = 2;
            str += 2;
        }
        else
        {
            base = 8;
            ++str;
        }
    }

    // Process each character according to the detected base.
    for (; *str; ++str)
    {
        char const c = *str;
        int digit = 0;

        if (c == '\'')
            continue;
        if (c >= '0' && c <= '9')
            digit = c - '0';
        else if (base == 16 && c >= 'a' && c <= 'f')
            digit = c - 'a' + 10;
        else if (base == 16 && c >= 'A' && c <= 'F')
            digit = c - 'A' + 10;
        else
            throw std::invalid_argument("Invalid character");

        if (digit >= base)
            throw std::invalid_argument("Digit out of range for the given base");

        res = res * base + digit;
    }
    return res;
}

inline namespace literals
{
consteval i8 operator""_i8(unsigned long long n) { return n; }
consteval i16 operator""_i16(unsigned long long n) { return n; }
consteval i32 operator""_i32(unsigned long long n) { return n; }
consteval i64 operator""_i64(unsigned long long n) { return n; }
consteval i128 operator""_i128(const char *str) { return fromString<i128>(str); }

consteval u8 operator""_u8(unsigned long long n) { return n; }
consteval u16 operator""_u16(unsigned long long n) { return n; }
consteval u32 operator""_u32(unsigned long long n) { return n; }
consteval u64 operator""_u64(unsigned long long n) { return n; }
consteval u128 operator""_u128(const char *str) { return fromString<u128>(str); }

inline cpp_int operator""_cppi(const char *str) { return cpp_int{str}; }
inline mpz_int operator""_Z(const char *str) { return mpz_int{str}; }
inline mpq_rational operator""_Q(const char *str) { return mpq_rational{str}; }
inline mpf_float operator""_R(const char *str) { return mpf_float{str}; }
} // namespace literals
} // namespace euler
