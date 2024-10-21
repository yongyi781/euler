#pragma once

#include "concepts.hpp"
#include <boost/multiprecision/integer.hpp>

inline namespace euler
{
#ifdef __SIZEOF_INT128__
using int128_t = __int128;
using uint128_t = unsigned __int128;
#else
using int128_t = boost::multiprecision::int128_t;
using uint128_t = boost::multiprecision::uint128_t;
#endif
using boost::multiprecision::int1024_t;
using boost::multiprecision::int256_t;
using boost::multiprecision::int512_t;
using int2048_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    2048, 2048, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using int4096_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    4096, 4096, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using int8192_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    8192, 8192, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using boost::container::static_vector;
using boost::multiprecision::cpp_int;
using boost::multiprecision::mpf_float;
using boost::multiprecision::mpq_rational;
using boost::multiprecision::mpz_int;
template <typename T> using PrimePower = std::pair<T, int>;

// template <typename T>
// struct FactorizationType
// {
//     using type = std::vector<PrimePower<T>>;
// };

// template <>
// struct FactorizationType<int>
// {
//     using type = static_vector<PrimePower<int>, 9>;
// };

// template <>
// struct FactorizationType<unsigned int>
// {
//     using type = static_vector<PrimePower<unsigned int>, 9>;
// };

// template <>
// struct FactorizationType<int64_t>
// {
//     using type = static_vector<PrimePower<int64_t>, 15>;
// };

// template <>
// struct FactorizationType<uint64_t>
// {
//     using type = static_vector<PrimePower<uint64_t>, 15>;
// };

// template <>
// struct FactorizationType<int128_t>
// {
//     using type = static_vector<PrimePower<int128_t>, 25>;
// };

// template <>
// struct FactorizationType<int256_t>
// {
//     using type = static_vector<PrimePower<int256_t>, 43>;
// };

// template <>
// struct FactorizationType<int512_t>
// {
//     using type = static_vector<PrimePower<int512_t>, 75>;
// };

template <typename T> using Factorization = std::vector<PrimePower<T>>;

/// Returns a base raised to an integer power. The type of the base needs a multiplication operation
/// defined on it.
template <typename T, integral2 U, std::invocable<T, T> BinaryOp>
constexpr T pow(const T &base, U exponent, const T &identity, BinaryOp op)
{
    if (exponent == 0 || base == identity)
        return identity;
    if (exponent == 1)
        return base;
    if constexpr ((integral2<T> || std::floating_point<T>) && std::is_same_v<BinaryOp, std::multiplies<>>)
        if (base == -identity)
            return exponent % 2 == 0 ? identity : -identity;
    assert(exponent >= 0);

    T x = identity;
    T y = base;
    while (true)
    {
        if (exponent & 1)
        {
            x = op(x, y);
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
        base *base;
        T(1);
    }
{
    return pow(base, exponent, T(1), std::multiplies{});
}

/// @brief Computes the correct modulo operation in the case of negative values.
/// @param a A number.
/// @param modulus A number.
/// @return A number between 0 and b - 1.
template <integral2 T, integral2 Tm>
    requires(!boost::multiprecision::is_number_expression<T>::value &&
             !boost::multiprecision::is_number_expression<Tm>::value)
constexpr std::common_type_t<T, Tm> mod(const T &a, Tm modulus)
{
    using Tp = std::common_type_t<T, Tm>;
    if (a >= 0)
        return Tp(a) < Tp(modulus) ? Tp(a) : Tp(a % modulus);
    return Tp((a + 1) % modulus + modulus - 1);
}

// template <typename T, size_t M, size_t N>
// constexpr Matrix<T, M, N> mod(const Matrix<T, M, N> &m, integral2 auto modulus);

/// Checked safe version of modmul.
template <integral2 Ta, integral2 Tb, integral2 Tm>
constexpr std::common_type_t<Ta, Tb, Tm> modmul(const Ta &a, const Tb &b, const Tm &m)
{
    using T = std::common_type_t<Ta, Tb, Tm>;
    if constexpr (std::numeric_limits<T>::digits == std::numeric_limits<int>::max())
    {
        // We're in the cpp_int or mpz_int case.
        return T(a * b % m);
    }
    else
    {
        using Td = boost::multiprecision::detail::double_integer<T>::type;
        if constexpr (requires(Td result) { __builtin_mul_overflow(a, b, &result); })
        {
            Td result{};
            __builtin_mul_overflow(a, b, &result);
            return T(result % m);
        }
        else
        {
            return T(Td(a) * Td(b) % Td(m));
        }
    }
}

/// Returns whether `a * b <= c`, and always returns false if the multiplication overflows.
template <integral2 T, integral2 U, integral2 V> constexpr bool mulLeq(T a, U b, V c)
{
    V x{};
#if defined(__clang__) || defined(__GNUC__)
    return !__builtin_mul_overflow(a, b, &x) && x <= c;
#else
    return !std::_Mul_overflow(V(a), V(b), x) && x <= c;
#endif
}

/// Computes the integral square root of a number.
template <integral2 T> constexpr T isqrt(T n)
{
    if (!std::integral<T> || std::is_constant_evaluated())
        return boost::multiprecision::sqrt(n);
    T x = (T)sqrt(n);
    return x * x > n ? x - 1 : x;
}

/// Returns the square root of an integer, if it is a square. Otherwise, returns none.
template <integral2 T> constexpr std::optional<T> sqrtIfSquare(T n)
{
    auto s = isqrt(n);
    return s * s == n ? std::optional{s} : std::nullopt;
}

/// Returns the square root of a rational number, if it is a square. Otherwise, returns none.
template <integral2 T> constexpr std::optional<boost::rational<T>> sqrt(const boost::rational<T> &r)
{
    T num = isqrt(r.numerator());
    T denom = isqrt(r.denominator());
    if (num * num != r.numerator() || denom * denom != r.denominator())
        return std::nullopt;
    return {{num, denom}};
}

template <integral2 T> constexpr bool isSquare(const T &n)
{
    T a = isqrt(n);
    return a * a == n;
}

template <integral2 T> constexpr bool isSquare(const boost::rational<T> &r)
{
    return isSquare(r.numerator()) && isSquare(r.denominator());
}

/// Calculates ⌊a / b⌋, and works for negative values too.
template <integral2 Ta, integral2 Tb> constexpr Ta floorDiv(const Ta &a, const Tb &b)
{
    auto d = a / b;
    return d * b == a ? Ta(d) : Ta(d - ((a < 0) ^ (b < 0)));
}

/// Calculates ⌈a / b⌉, and works for negative values too.
template <integral2 Ta, integral2 Tb> constexpr Ta ceilDiv(const Ta &a, const Tb &b)
{
    return floorDiv<Ta, Tb>(a + b - 1, b);
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

/// Greatest common divisor of a range.
template <std::ranges::range Range> std::ranges::range_value_t<Range> gcd(Range &&r)
{
    using T = std::ranges::range_value_t<Range>;
    using std::gcd;
    T g = 0;
    for (auto &&x : r)
    {
        g = gcd(g, x);
        if (g == 1)
            break;
    }
    return g;
}

/// Greatest common divisor, specialized to an initializer list.
template <integral2 T> T gcd(std::initializer_list<T> l)
{
    using std::gcd;
    T g = 0;
    for (auto &&x : l)
    {
        g = gcd(g, x);
        if (g == 1)
            break;
    }
    return g;
}

/// Least common multiple of a range.
template <std::ranges::range Range> std::ranges::range_value_t<Range> lcm(Range &&r)
{
    using T = std::ranges::range_value_t<Range>;
    using std::lcm;
    T l = 1;
    for (auto &&x : r)
        l = lcm(l, x);
    return l;
}

/// Fast binomial coefficient for small values. Guaranteed correct for `T = int64_t` only for `n` up to 60, and for
/// `T = uint64_t` only for `n` up to 62.
constexpr uint64_t binom(int n, int k)
{
    using std::min;
    if (k < 0 || k > n)
        return 0;
    k = min(k, n - k);
    uint64_t result = 1;
    for (int i = 0; i < k; ++i)
        result = result * (n - i) / (i + 1);
    return result;
}

// Function to check if a character is the start of a UTF-8 sequence
constexpr bool isUtf8Start(char c) { return (c & 0xC0) != 0x80; }

// Function to calculate the display width of a UTF-8 character
constexpr int charWidth(char32_t c)
{
    if (c == 0)
        return 0;
    if (c < 32 || (c >= 0x7F && c < 0xA0))
        return 0; // Control characters
    if (c >= 0x1100 &&
        (c <= 0x115F ||                                                               // Hangul Jamo
         c == 0x2329 || c == 0x232A || (c >= 0x2E80 && c <= 0xA4CF && c != 0x303F) || // CJK ... Yi
         (c >= 0xAC00 && c <= 0xD7A3) ||                                              // Hangul Syllables
         (c >= 0xF900 && c <= 0xFAFF) ||                                              // CJK Compatibility Ideographs
         (c >= 0xFE10 && c <= 0xFE19) ||                                              // Vertical forms
         (c >= 0xFE30 && c <= 0xFE6F) ||                                              // CJK Compatibility Forms
         (c >= 0xFF00 && c <= 0xFF60) ||                                              // Fullwidth Forms
         (c >= 0xFFE0 && c <= 0xFFE6) || (c >= 0x20000 && c <= 0x2FFFD) || (c >= 0x30000 && c <= 0x3FFFD)))
    {
        return 2;
    }
    return 1;
}

// Function to calculate display width of a string, ignoring ANSI escape codes
constexpr size_t displayWidth(std::string_view str)
{
    size_t width = 0;
    bool in_escape = false;
    const auto *it = str.begin();

    while (it != str.end())
    {
        if (*it == '\033')
        {
            in_escape = true;
            ++it;
            continue;
        }

        if (in_escape)
        {
            if (std::isalpha(static_cast<uint8_t>(*it)))
                in_escape = false;
            ++it;
            continue;
        }

        if ((*it & 0x80) == 0)
        {
            width += charWidth(*it);
            ++it;
        }
        else
        {
            char32_t cp = 0;
            int bytes = 0;
            if ((*it & 0xE0) == 0xC0)
            {
                cp = (*it & 0x1F);
                bytes = 2;
            }
            else if ((*it & 0xF0) == 0xE0)
            {
                cp = (*it & 0x0F);
                bytes = 3;
            }
            else if ((*it & 0xF8) == 0xF0)
            {
                cp = (*it & 0x07);
                bytes = 4;
            }
            else
            {
                // Invalid UTF-8, skip
                ++it;
                continue;
            }

            for (int i = 1; i < bytes && it + i != str.end(); ++i)
                cp = (cp << 6) | (*(it + i) & 0x3F);

            width += charWidth(cp);
            std::advance(it, bytes);
        }
    }

    return width;
}

/// Measures the number of characters that would be output by `o << x`.
template <typename T, typename CharT = char, typename Traits = std::char_traits<CharT>>
std::string toStringWithFlags(const T &x, const std::basic_ostream<CharT, Traits> &o = std::cout);

// Forward declarations of printing operators.

#ifdef __SIZEOF_INT128__
template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const int128_t &x);
template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const uint128_t &x);
template <typename CharT, typename Traits>
std::basic_istream<CharT, Traits> &operator>>(std::basic_istream<CharT, Traits> &is, int128_t &x);
template <typename CharT, typename Traits>
std::basic_istream<CharT, Traits> &operator>>(std::basic_istream<CharT, Traits> &is, uint128_t &x);
#endif

template <typename CharT, typename Traits, typename T, typename U>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const std::pair<T, U> &v);

template <typename CharT, typename Traits, typename... Args>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, std::tuple<Args...> const &t);

template <typename CharT, typename Traits, typename T>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const std::optional<T> &x);

template <typename CharT, typename Traits, typename T>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, std::nullopt_t /*unused*/);

template <typename CharT, typename Traits, std::ranges::range Range>
    requires(!is_string<Range>)
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Range &r);

template <integral2 TMod> class mod_plus
{
  public:
    explicit constexpr mod_plus(TMod modulus) : modulus(std::move(modulus)) {}

    template <typename T, typename U>
    constexpr decltype(auto) operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) + std::forward<U>(b)))
    {
        return (std::forward<T>(a) + std::forward<U>(b)) % modulus;
    }

  private:
    TMod modulus;
};

/** Check for overflow potential! */
template <typename TMod> class mod_multiplies
{
  public:
    explicit constexpr mod_multiplies(TMod modulus) : modulus(std::move(modulus)) {}

    template <typename T, typename U>
    constexpr decltype(auto) operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) * std::forward<U>(b)))
    {
        return std::forward<T>(a) * std::forward<U>(b) % modulus;
    }

  private:
    TMod modulus;
};

/** Check for overflow potential! */
template <integral2 TMod> class mod_multiplies_safe
{
  public:
    explicit constexpr mod_multiplies_safe(TMod modulus) : modulus(std::move(modulus)) {}

    template <typename T, typename U>
    constexpr std::common_type_t<T, U, TMod> operator()(T &&a, U &&b) const
        noexcept(noexcept(std::forward<T>(a) * std::forward<U>(b)))
    {
        return modmul(std::forward<T>(a), std::forward<U>(b), modulus);
    }

  private:
    TMod modulus;
};

struct minimum
{
    template <typename T, typename U>
    constexpr auto operator()(T &&a, U &&b) const
        noexcept(noexcept(std::forward<T>(a) + std::forward<U>(b))) -> decltype(std::min(std::forward<T>(a),
                                                                                         std::forward<U>(b)))
    {
        return std::min(std::forward<T>(a), std::forward<U>(b));
    }
};

struct maximum
{
    template <typename T, typename U>
    constexpr auto operator()(T &&a, U &&b) const
        noexcept(noexcept(std::forward<T>(a) + std::forward<U>(b))) -> decltype(std::max(std::forward<T>(a),
                                                                                         std::forward<U>(b)))
    {
        return std::max(std::forward<T>(a), std::forward<U>(b));
    }
};

// Removes the dependency on boost::counting_iterator.
template <typename T> struct counting_iterator
{
    using iterator_category = std::random_access_iterator_tag;
    using iterator_concept = std::random_access_iterator_tag;
    using difference_type = std::iter_difference_t<T>;
    using pointer = const T *;
    using reference = const T &;
    using value_type = T;

    counting_iterator() = default;
    constexpr counting_iterator(T p) : _value(std::move(p)) {}

    constexpr T operator*() const { return _value; }

    constexpr counting_iterator &operator++()
    {
        ++_value;
        return *this;
    }

    constexpr counting_iterator operator++(int)
    {
        counting_iterator tmp = *this;
        ++(*this);
        return tmp;
    }

    constexpr counting_iterator &operator+=(T i)
    {
        _value += i;
        return *this;
    }

    constexpr counting_iterator operator+(const difference_type other) const { return {T(_value + other)}; }
    friend constexpr counting_iterator operator+(const difference_type value, const counting_iterator &other)
    {
        return {T(other + value)};
    }

    constexpr counting_iterator &operator--()
    {
        --_value;
        return *this;
    }
    constexpr counting_iterator operator--(int)
    {
        counting_iterator tmp = *this;
        --(*this);
        return tmp;
    }
    constexpr counting_iterator &operator-=(T i)
    {
        _value -= i;
        return *this;
    }
    constexpr difference_type operator-(const counting_iterator &other) const { return _value - other._value; }

    constexpr counting_iterator operator-(const difference_type other) const { return _value - other; }
    friend constexpr counting_iterator operator-(const difference_type value, const counting_iterator &other)
    {
        return other - value;
    }

    constexpr T operator[](difference_type i) const { return _value + i; }

    constexpr std::strong_ordering operator<=>(const counting_iterator &other) const
    {
        if (_value == other._value)
            return std::strong_ordering::equal;
        return _value < other._value ? std::strong_ordering::less : std::strong_ordering::greater;
    }

    constexpr bool operator==(const counting_iterator &other) const { return _value == other._value; }
    constexpr bool operator!=(const counting_iterator &other) const { return _value != other._value; }

  private:
    T _value;
};

#ifdef __cpp_multidimensional_subscript
/// 2D vector.
template <typename T> class vector2d
{
  public:
    using value_type = T;
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

    vector2d() = default;
    constexpr vector2d(size_t rows, size_t columns) : _rows(rows), _columns(columns), _data(rows * columns) {}
    constexpr vector2d(size_t rows, size_t columns, const T &value) : vector2d(rows, columns)
    {
        std::ranges::fill(_data, value);
    }

    /// Gets the linear index corresponding to a coordinate `(i, j)`.
    [[nodiscard]] constexpr int toIndex(int i, int j) const noexcept { return _columns * i + j; }

    /// `i`th row, `j`th column.
    constexpr const T &operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr T &operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

    [[nodiscard]] constexpr size_t rows() const noexcept { return _rows; }
    [[nodiscard]] constexpr size_t columns() const noexcept { return _columns; }

    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }

    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    constexpr void resize(size_t rows, size_t columns)
    {
        _rows = rows;
        _columns = columns;
        _data.resize(columns * rows);
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const vector2d<T> &v)
    {
        int maxWidth = o.width();
        for (auto &&x : v._data)
            maxWidth = std::max(maxWidth, (int)toStringWithFlags(x, o).size());
        std::basic_ostringstream<CharT, Traits> ss;
        ss.flags(o.flags());
        ss.imbue(o.getloc());
        ss.precision(o.precision());
        ss << v._columns << "×" << v._rows << " vector2d:\n";
        for (size_t i = 0; i < v.rows(); ++i)
        {
            for (size_t j = 0; j < v.columns(); ++j)
                ss << std::setw(maxWidth + 1) << v[i, j];
            ss << '\n';
        }
        return o << std::move(ss).str();
    }

  private:
    size_t _rows = 0;
    size_t _columns = 0;
    std::vector<T> _data;
};

// Triangular vector. Consists of indices `(i, j)` where `j ≤ i`.
template <typename T> class triangular_vector
{
  public:
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

    triangular_vector() = default;
    constexpr triangular_vector(size_t height) : _height(height), _data(toIndex(height, 0)) {}

    /// `i`th row, `j`th column.
    constexpr const T &operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr T &operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

    [[nodiscard]] constexpr size_t height() const noexcept { return _height; }

    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }

    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    constexpr void resize(size_t newHeight)
    {
        _height = newHeight;
        _data.resize(toIndex(newHeight, 0));
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const triangular_vector &v)
    {
        for (size_t i = 0; i < v.height(); ++i)
        {
            for (size_t j = 0; j < v.width(); ++j)
            {
                if (j != 0)
                    o << ", ";
                o << v[i, j];
            }
            o << '\n';
        }
        return o;
    }

  private:
    size_t _height = 0;
    std::vector<T> _data;

    static constexpr size_t toIndex(size_t i, size_t j) noexcept { return (i * (i + 1) / 2) + j; }
};

// Triangular array. Consists of indices `(i, j)` where `j ≤ i`.
template <typename T, size_t N> class triangular_array
{
  public:
    using iterator = std::array<T, N *(N + 1) / 2>::iterator;
    using const_iterator = std::array<T, N *(N + 1) / 2>::const_iterator;

    /// `i`th row, `j`th column.
    constexpr const T &operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr T &operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }
    [[nodiscard]] constexpr size_t height() const noexcept { return N; }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }

    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const triangular_array &v)
    {
        for (size_t i = 0; i < v.height(); ++i)
        {
            for (size_t j = 0; j < v.width(); ++j)
            {
                if (j != 0)
                    o << ", ";
                o << v[i, j];
            }
            o << '\n';
        }
        return o;
    }

  private:
    std::array<T, N *(N + 1) / 2> _data;

    static constexpr size_t toIndex(size_t i, size_t j) noexcept
    {
        return j <= i ? (i * (i + 1) / 2) + j : (j * (j + 1) / 2) + i;
    }
};
#endif

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

#ifdef DEBUG
// String routines for natvis
std::string to_string_cpp_int(const cpp_int &x) { return x.backend().size() == 0 ? "uninitialized" : x.str(); }
std::string to_string_mpz_int(const mpz_int &x) { return x.str(); }
#endif
} // namespace euler

namespace std
{
template <> inline int128_t gcd(int128_t __m, int128_t __n) noexcept { return boost::integer::gcd(__m, __n); }
template <> inline int128_t lcm(int128_t __m, int128_t __n) noexcept { return boost::integer::lcm(__m, __n); }
template <> inline uint128_t gcd(uint128_t __m, uint128_t __n) noexcept { return boost::integer::gcd(__m, __n); }
template <> inline uint128_t lcm(uint128_t __m, uint128_t __n) noexcept { return boost::integer::lcm(__m, __n); }
} // namespace std
