#pragma once

#include "matrix.hpp"
#include <ostream>

namespace euler
{
/// Represents an affine linear transformation a * x + b.
template <typename T> struct Aff
{
    T a, b;

    Aff(T a = 1, T b = 0) : a(a), b(b) {}

    [[nodiscard]] constexpr T operator()(T x) const { return a * x + b; }

    constexpr Aff &operator+=(const Aff &o)
    {
        a += o.a;
        b += o.b;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator+(Aff left, const Aff &right)
    {
        left += right;
        return left;
    }

    constexpr Aff &operator-=(const Aff &o)
    {
        a -= o.a;
        b -= o.b;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator-(Aff left, const Aff &right)
    {
        left -= right;
        return left;
    }

    constexpr Aff &operator*=(const Aff &o)
    {
        b += a * o.b;
        a *= o.a;
        return *this;
    }
    [[nodiscard]] constexpr friend Aff operator*(Aff left, const Aff &right)
    {
        left *= right;
        return left;
    }

    /// Returns the inverse of this affine linear transformation.
    [[nodiscard]] constexpr Aff operator~() const { return inverse(); }

    /// Returns the inverse of this affine linear transformation.
    [[nodiscard]] constexpr Aff inverse() const
    {
        T const inv_a = T(1) / a;
        return Aff(inv_a, -b * inv_a);
    }

    /// Returns the matrix representation of this affine linear transformation.
    [[nodiscard]] constexpr Matrix<T, 2> mat() const { return Matrix<T, 2>{{a, b}, {1, 0}}; }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &os, const Aff &a)
    {
        return os << a.a << "*x + " << a.b;
    }
};
} // namespace euler
