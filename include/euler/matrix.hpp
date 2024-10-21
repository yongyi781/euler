#pragma once

#include "decls.hpp"
#include "modular_arithmetic.hpp"
#include <array>

inline namespace euler
{
template <typename T, size_t M, size_t N = M> class Matrix;

/// A stack-allocated vector class.
template <typename T, size_t N> class Vector
{
  public:
    using value_type = T;
    using iterator = std::array<T, N>::iterator;
    using const_iterator = std::array<T, N>::const_iterator;
    static constexpr size_t length = N;

    [[nodiscard]] Vector() = default;
    [[nodiscard]] constexpr Vector(std::array<T, N> data) : _data(std::move(data)) {}
    template <typename... U> [[nodiscard]] constexpr Vector(U... data) : _data{data...} {}

    [[nodiscard]] constexpr static Vector zero() { return {}; }

    [[nodiscard]] constexpr T &operator[](size_t i) { return _data[i]; }
    [[nodiscard]] constexpr const T &operator[](size_t i) const { return _data[i]; }

    constexpr Vector &operator+=(const Vector &other)
    {
        for (size_t i = 0; i < N; ++i)
            _data[i] += other._data[i];
        return *this;
    }
    [[nodiscard]] constexpr Vector operator+(const Vector &other) const { return Vector{*this} += other; }

    constexpr Vector &operator-=(const Vector &other)
    {
        for (size_t i = 0; i < N; ++i)
            _data[i] -= other._data[i];
        return *this;
    }
    [[nodiscard]] constexpr Vector operator-(const Vector &other) const { return Vector{*this} -= other; }

    constexpr Vector &operator*=(T t)
    {
        for (auto &&x : _data)
            x *= t;
        return *this;
    }
    [[nodiscard]] constexpr Vector operator*(T t) const { return Vector{*this} *= t; }
    [[nodiscard]] friend constexpr Vector operator*(T t, const Vector &m) { return m * t; }

    constexpr Vector &operator/=(T t)
    {
        for (auto &&x : _data)
            x /= t;
        return *this;
    }
    [[nodiscard]] constexpr Vector operator/(T t) { return Vector{*this} /= t; }

    constexpr Vector &operator%=(T t)
    {
        for (auto &&x : _data)
            x %= t;
        return *this;
    }
    [[nodiscard]] constexpr Vector operator%(T t) const { return Vector{*this} %= t; }

    [[nodiscard]] constexpr Vector operator-() const
    {
        Vector result{};
        for (size_t i = 0; i < N; ++i)
            result._data[i] = -_data[i];
        return result;
    }

    [[nodiscard]] constexpr Vector operator+() const { return *this; }

    constexpr bool operator==(const Vector &other) const { return _data == other._data; }
    constexpr bool operator!=(const Vector &other) const { return _data != other._data; }

    constexpr auto operator<=>(const Vector &other) const { return _data <=> other._data; }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Vector &x)
    {
        return o << x._data;
    }

    [[nodiscard]] constexpr std::array<T, N> &data() { return _data; }
    [[nodiscard]] constexpr const std::array<T, N> &data() const { return _data; }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }
    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    /// Makes a row vector, i.e. a 1 × N matrix.
    [[nodiscard]] constexpr Matrix<T, 1, N> transpose() const
    {
        Matrix<T, 1, N> result{};
        for (size_t i = 0; i < N; ++i)
            result[0, i] = _data[i];
        return result;
    }

    /// Makes a diagonal matrix out of this vectdor.
    [[nodiscard]] constexpr Matrix<T, N, N> toDiagonal() const
    {
        Matrix<T, N, N> result{};
        for (size_t i = 0; i < N; ++i)
            result[i, i] = _data[i];
        return result;
    }

    /// Dot product.
    [[nodiscard]] constexpr T dot(const Vector &other) const
    {
        T result{};
        for (size_t i = 0; i < N; ++i)
            result += _data[i] * other._data[i];
        return result;
    }

    // 2D cross product, outputs a scalar.
    [[nodiscard]] constexpr T cross(const Vector &other) const
        requires(N == 2)
    {
        auto &&[a1, a2] = _data;
        auto &&[b1, b2] = other._data;
        return a1 * b2 - a2 * b1;
    }

    // 3D cross product.
    [[nodiscard]] constexpr Vector cross(const Vector &other) const
        requires(N == 3)
    {
        auto &&[a1, a2, a3] = _data;
        auto &&[b1, b2, b3] = other._data;
        return {a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1};
    }

  private:
    std::array<T, N> _data;
};

template <typename T, typename... U>
    requires(std::is_same_v<T, U> && ...)
Vector(T, U...) -> Vector<T, 1 + sizeof...(U)>;
template <std::size_t I, typename T, std::size_t N>
    requires(I < N)
[[nodiscard]] constexpr T &get(Vector<T, N> &v) noexcept
{
    return v[I];
}

template <std::size_t I, typename T, std::size_t N>
    requires(I < N)
[[nodiscard]] constexpr T &&get(Vector<T, N> &&v) noexcept
{
    return std::move(std::get<I>(v));
}

template <std::size_t I, typename T, std::size_t N>
    requires(I < N)
[[nodiscard]] constexpr const T &get(const Vector<T, N> &v) noexcept
{
    return v[I];
}

template <std::size_t I, typename T, std::size_t N>
    requires(I < N)
[[nodiscard]] constexpr const T &&get(const Vector<T, N> &&v) noexcept
{
    return std::move(std::get<I>(v));
}

/// A stack-allocated matrix with compile-time constant size.
template <typename T, size_t M, size_t N> class Matrix
{
  public:
    [[nodiscard]] Matrix() = default;
    [[nodiscard]] constexpr Matrix(std::array<std::array<T, N>, M> data) : _data(std::move(data)) {}
    template <std::invocable<size_t, size_t> Fun> [[nodiscard]] constexpr Matrix(Fun f)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                _data[i][j] = f(i, j);
    }

    /// Makes a multiple of the identity matrix.
    template <integral2 U>
        requires(M == N && std::convertible_to<U, T>)
    [[nodiscard]] constexpr Matrix(U scalar)
    {
        for (size_t i = 0; i < N; ++i)
            _data[i][i] = scalar;
    }

    [[nodiscard]] constexpr static Matrix zero() { return {}; }

    [[nodiscard]] constexpr static Matrix identity()
        requires(M == N)
    {
        return {T(1)};
    }

    /// Elementary matrix with 1 in the (i, j) entry.
    [[nodiscard]] constexpr static Matrix elementary(size_t i, size_t j)
        requires(M == N)
    {
        assert(i >= 0 && i < M && j >= 0 && j < N);
        auto result = identity();
        result._data[i][j] = 1;
        return result;
    }

#ifdef __cpp_multidimensional_subscript
    /// Accesses the (i, j) entry, by reference.
    [[nodiscard]] constexpr T &operator[](size_t i, size_t j) { return _data[i][j]; }
    /// Accesses the (i, j) entry, by const reference.
    [[nodiscard]] constexpr const T &operator[](size_t i, size_t j) const { return _data[i][j]; }
#endif

    /// Accesses the ith row, by reference.
    [[nodiscard]] constexpr std::array<T, N> &operator[](size_t i) { return _data[i]; }
    /// Accesses the ith row, by const reference.
    [[nodiscard]] constexpr const std::array<T, N> &operator[](size_t i) const { return _data[i]; }

    constexpr Matrix &operator+=(const Matrix &other)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                _data[i][j] += other._data[i][j];
        return *this;
    }
    [[nodiscard]] constexpr Matrix operator+(const Matrix &other) const { return Matrix{*this} += other; }

    constexpr Matrix &operator-=(const Matrix &other)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                _data[i][j] -= other._data[i][j];
        return *this;
    }
    [[nodiscard]] constexpr Matrix operator-(const Matrix &other) const { return Matrix{*this} -= other; }

    template <size_t P> [[nodiscard]] constexpr Matrix<T, M, P> operator*(const Matrix<T, N, P> &other) const
    {
        Matrix<T, M, P> result{};
        for (size_t i = 0; i < M; ++i)
            for (size_t k = 0; k < N; ++k)
                for (size_t j = 0; j < P; ++j)
                    result[i, j] += _data[i][k] * other[k, j];
        return result;
    }
    template <size_t P> constexpr Matrix<T, M, P> operator*=(const Matrix<T, N, P> &other)
    {
        return *this = *this * other;
    }

    [[nodiscard]] constexpr Vector<T, M> operator*(const Vector<T, N> &v) const
    {
        Vector<T, M> result{};
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result[i] += _data[i][j] * v[j];
        return result;
    }

    constexpr Matrix &operator*=(T t)
    {
        for (auto &&row : _data)
            for (auto &&x : row)
                x *= t;
        return *this;
    }
    [[nodiscard]] constexpr Matrix operator*(T t) const { return Matrix{*this} *= t; }
    [[nodiscard]] friend constexpr Matrix operator*(T t, const Matrix &m) { return m * t; }

    constexpr Matrix &operator/=(T t)
    {
        for (auto &&row : _data)
            for (auto &&x : row)
                x /= t;
        return *this;
    }
    [[nodiscard]] constexpr Matrix operator/(T t) { return Matrix{*this} /= t; }

    constexpr Matrix &operator%=(T t)
    {
        for (auto &&row : _data)
            for (auto &&x : row)
                x %= t;
        return *this;
    }
    [[nodiscard]] constexpr Matrix operator%(T t) const { return Matrix{*this} %= t; }

    [[nodiscard]] constexpr Matrix operator-() const
    {
        Matrix result{};
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result._data[i][j] = -_data[i][j];
        return result;
    }

    [[nodiscard]] constexpr Matrix operator+() const { return *this; }

    constexpr bool operator==(const Matrix &other) const { return _data == other._data; }
    constexpr bool operator!=(const Matrix &other) const { return _data != other._data; }

    constexpr std::strong_ordering operator<=>(const Matrix &other) const { return _data <=> other._data; }

    [[nodiscard]] std::array<std::array<T, N>, M> data() const { return _data; }

    [[nodiscard]] constexpr Matrix<T, N, M> transpose() const
    {
        Matrix<T, N, M> result{};
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result[j, i] = _data[i][j];
        return result;
    }

    /// Returns the power of an matrix to an integer.
    [[nodiscard]] constexpr Matrix pow(integral2 auto exponent) const
        requires(M == N)
    {
        return euler::pow(*this, std::move(exponent), identity(), std::multiplies{});
    }

    /// Returns the power of an matrix to an integer.
    [[nodiscard]] constexpr Matrix powm(integral2 auto exponent, integral2 auto modulus)
        requires(M == N)
    {
        return euler::powm(*this, std::move(exponent), std::move(modulus), identity());
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Matrix &m)
    {
        int maxWidth = o.width();
        for (auto &&xs : m._data)
            for (auto &&x : xs)
                maxWidth = std::max(maxWidth, (int)toStringWithFlags(x, o).size());
        std::basic_ostringstream<CharT, Traits> ss;
        ss.flags(o.flags());
        ss.imbue(o.getloc());
        ss.precision(o.precision());
        ss << M << "×" << N << " matrix:\n";
        for (auto &&xs : m._data)
        {
            for (auto &&x : xs)
                ss << std::setw(maxWidth + 1) << x;
            ss << '\n';
        }
        return o << std::move(ss).str();
    }

  private:
    std::array<std::array<T, N>, M> _data{};
};

/// A stack-allocated symmetric matrix with compile-time constant size.
template <typename T, size_t N> class SymmetricMatrix
{
  public:
    static constexpr size_t size = N * (N + 1) / 2;

    [[nodiscard]] SymmetricMatrix() = default;
    [[nodiscard]] constexpr SymmetricMatrix(std::array<T, size> data) : _data(std::move(data)) {}

    /// Makes a multiple of the identity symmetric matrix.
    [[nodiscard]] constexpr SymmetricMatrix(T scalar) : _data{}
    {
        for (size_t i = 0; i < N; ++i)
            (*this)[i, i] = scalar;
    }

    [[nodiscard]] constexpr static SymmetricMatrix zero() { return {}; }

    [[nodiscard]] constexpr static SymmetricMatrix identity() { return {T(1)}; }

    /// Elementary symmetric matrix with 1 in the (i, j) and (j, i) entry.
    [[nodiscard]] constexpr static SymmetricMatrix elementary(size_t i, size_t j)
    {
        auto result = identity();
        result[i, j] = 1;
        return result;
    }

#ifdef __cpp_multidimensional_subscript

    /// Accesses the (i, j) entry, by reference.
    [[nodiscard]] constexpr T &operator[](size_t i, size_t j)
    {
        return j <= i ? _data[i * (i + 1) / 2 + j] : _data[j * (j + 1) / 2 + i];
    }

    /// Accesses the (i, j) entry, by const reference.
    [[nodiscard]] constexpr const T &operator[](size_t i, size_t j) const
    {
        return j <= i ? _data[i * (i + 1) / 2 + j] : _data[j * (j + 1) / 2 + i];
    }
#endif

    constexpr SymmetricMatrix &operator+=(const SymmetricMatrix &other)
    {
        for (size_t i = 0; i < size; ++i)
            _data[i] += other._data[i];
        return *this;
    }
    [[nodiscard]] constexpr SymmetricMatrix operator+(const SymmetricMatrix &other) const
    {
        return SymmetricMatrix{*this} += other;
    }

    constexpr SymmetricMatrix &operator-=(const SymmetricMatrix &other)
    {
        for (size_t i = 0; i < size; ++i)
            _data[i] += other._data[i];
        return *this;
    }
    [[nodiscard]] constexpr SymmetricMatrix operator-(const SymmetricMatrix &other) const
    {
        return SymmetricMatrix{*this} -= other;
    }

    [[nodiscard]] constexpr SymmetricMatrix operator*(const SymmetricMatrix &other) const
    {
        SymmetricMatrix result{};
        for (size_t i = 0; i < N; ++i)
            for (size_t k = 0; k < N; ++k)
                for (size_t j = 0; j <= i; ++j)
                    result[i, j] += (*this)[i, k] * other[k, j];
        return result;
    }
    constexpr SymmetricMatrix operator*=(const SymmetricMatrix &other) { return *this = *this * other; }

    [[nodiscard]] constexpr Vector<T, N> operator*(const Vector<T, N> &v) const
    {
        Vector<T, N> result{};
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                result[i] += _data[i, j] * v[j];
        return result;
    }

    constexpr SymmetricMatrix &operator*=(T t)
    {
        for (auto &&x : _data)
            x *= t;
        return *this;
    }
    [[nodiscard]] constexpr SymmetricMatrix operator*(T t) const { return SymmetricMatrix{*this} *= t; }
    [[nodiscard]] friend constexpr SymmetricMatrix operator*(T t, const SymmetricMatrix &m) { return m * t; }

    constexpr SymmetricMatrix &operator/=(T t)
    {
        for (auto &&x : _data)
            x /= t;
        return *this;
    }
    [[nodiscard]] constexpr SymmetricMatrix operator/(T t) { return SymmetricMatrix{*this} /= t; }

    constexpr SymmetricMatrix &operator%=(T t)
    {
        for (auto &&x : _data)
            x %= t;
        return *this;
    }
    [[nodiscard]] constexpr SymmetricMatrix operator%(T t) const { return SymmetricMatrix{*this} %= t; }

    [[nodiscard]] constexpr SymmetricMatrix operator-() const
    {
        SymmetricMatrix result{};
        for (auto &&x : result._data)
            x = -x;
        return result;
    }

    [[nodiscard]] constexpr SymmetricMatrix operator+() const { return *this; }

    constexpr bool operator==(const SymmetricMatrix &other) const { return _data == other._data; }
    constexpr bool operator!=(const SymmetricMatrix &other) const { return _data != other._data; }

    constexpr std::strong_ordering operator<=>(const SymmetricMatrix &other) const { return _data <=> other._data; }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const SymmetricMatrix &x)
    {
        return o << x._data;
    }

    [[nodiscard]] std::array<T, size> data() const { return _data; }

    [[nodiscard]] constexpr SymmetricMatrix transpose() const { return *this; }

    /// Returns the power of an SymmetricMatrix to an integer.
    [[nodiscard]] constexpr SymmetricMatrix pow(integral2 auto exponent) const
    {
        return ::pow(*this, std::move(exponent), identity(), std::multiplies{});
    }

    /// Returns the power of an SymmetricMatrix to an integer.
    [[nodiscard]] constexpr SymmetricMatrix powm(integral2 auto exponent, integral2 auto modulus)
    {
        return ::powm(*this, std::move(exponent), std::move(modulus), identity());
    }

  private:
    std::array<T, size> _data;
};

template <typename T> using Matrix2 = Matrix<T, 2>;
using Matrix2i = Matrix<int, 2>;
using Matrix2ll = Matrix<int64_t, 2>;
using Matrix2i128 = Matrix<int128_t, 2>;
template <typename T> using Matrix3 = Matrix<T, 3>;
using Matrix3i = Matrix<int, 3>;
using Matrix3ll = Matrix<int64_t, 3>;
using Matrix3i128 = Matrix<int128_t, 3>;

/// Returns the coefficient of `x^n` in `(∑ i ∈ [0..M-1], a[i] * x^i) / (∑ i ∈ [0..N-1], b[i] * x^i)`.
/// Requires `M < N`.
template <integral2 T, integral2 U, size_t M, size_t N>
    requires(M < N)
auto seriesCoefficient(const std::array<T, M> &a, const std::array<T, N> &b, U n)
{
    auto b0inv = T(1) / b[0];
    Matrix<T, N - 1> m{};
    for (size_t j = 0; j < N - 1; ++j)
        m[N - 2][j] = -b[N - j - 1] * b0inv;
    for (size_t i = 1; i < N - 1; ++i)
        m[i - 1][i] = 1;

    Vector<T, N - 1> v{};
    for (size_t i = 0; i < N - 1; ++i)
    {
        T t = i < a.size() ? a[i] : T(0);
        for (size_t j = 0; j < i; ++j)
            t = (t - b[i - j] * v[j]);
        v[i] = b0inv * t;
    }
    return (pow(m, n) * v)[0];
}

template <typename T, size_t N> constexpr size_t hash_value(const Vector<T, N> &v)
{
    return boost::hash_value(v.data());
}
} // namespace euler

namespace std
{
/// Partial specialization for Vector
template <typename T, size_t N> struct tuple_size<Vector<T, N>> : public integral_constant<size_t, N>
{
};

/// Partial specialization for Vector
template <size_t I, typename T, size_t N>
    requires(I < N)
struct tuple_element<I, Vector<T, N>>
{
    using type = T;
};
} // namespace std
