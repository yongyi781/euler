#pragma once

#include "matrix.hpp"
#include <ostream>

namespace euler
{
/// Represents an affine linear transformation (a * x + b) / c.
template <typename T> struct Aff
{
    T a, b, c;

    Aff(T a = 1, T b = 0, T c = 1) : a(a), b(b), c(c) {}

    [[nodiscard]] constexpr T operator()(T x) const { return (a * x + b) / c; }

    constexpr Aff &operator+=(const Aff &o)
    {
        a = a * o.c + o.a * c;
        b = b * o.c + o.b * c;
        c *= o.c;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator+(Aff left, const Aff &right)
    {
        left += right;
        return left;
    }

    constexpr Aff &operator-=(const Aff &o)
    {
        a = a * o.c - o.a * c;
        b = b * o.c - o.b * c;
        c *= o.c;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator-(Aff left, const Aff &right)
    {
        left -= right;
        return left;
    }

    constexpr Aff &operator*=(const Aff &o)
    {
        b = a * o.b + b * o.c;
        a *= o.a;
        c *= o.c;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator*(Aff left, const Aff &right)
    {
        left *= right;
        return left;
    }

    constexpr Aff &operator*=(T scalar)
    {
        a *= scalar;
        b *= scalar;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator*(Aff f, T scalar)
    {
        f *= scalar;
        return f;
    }
    [[nodiscard]] constexpr friend Aff operator*(T scalar, Aff f)
    {
        f *= scalar;
        return f;
    }

    /// Returns the inverse of this affine linear transformation.
    [[nodiscard]] constexpr Aff operator~() const { return inverse(); }

    [[nodiscard]] constexpr friend bool operator==(const Aff &left, const Aff &right)
    {
        return left.a * right.c == right.a * left.c && left.b * right.c == right.b * left.c;
    }

    /// Returns the inverse of this affine linear transformation.
    [[nodiscard]] constexpr Aff inverse() const { return {c, -b, a}; }

    /// Normalizes this affine linear transformation.
    [[nodiscard]] constexpr Aff normalize() const { return {a / c, b / c, 1}; }

    /// Returns the matrix representation of this affine linear transformation.
    [[nodiscard]] constexpr Matrix<T, 2> mat() const { return Matrix<T, 2>{{a, b}, {0, c}}; }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &os, const Aff &f)
    {
        if (f.c == 1)
            return os << f.a << "*x + " << f.b;
        return os << '(' << f.a << "*x + " << f.b << ")/" << f.c;
    }
};
} // namespace euler
