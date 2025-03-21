#pragma once

#include "concepts.hpp"
#include <boost/multiprecision/gmp.hpp>
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
using boost::multiprecision::uint1024_t;
using boost::multiprecision::uint256_t;
using boost::multiprecision::uint512_t;
using int2048_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    2048, 2048, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using int4096_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    4096, 4096, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using int8192_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    8192, 8192, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using boost::container::static_vector;
using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using boost::multiprecision::mpf_float;
using boost::multiprecision::mpq_rational;
using boost::multiprecision::mpz_int;

template <typename T> using PrimePower = std::pair<T, int>;

template <typename T> struct double_integer
{
    using type = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
        2 * std::numeric_limits<T>::digits, 2 * std::numeric_limits<T>::digits,
        boost::multiprecision::is_signed_number<T>::value ? boost::multiprecision::signed_magnitude
                                                          : boost::multiprecision::unsigned_magnitude>>;
};

template <std::integral T>
    requires(std::is_signed_v<T> && 2 * sizeof(T) <= sizeof(intmax_t))
struct double_integer<T>
{
    using type = boost::int_t<CHAR_BIT * 2 * sizeof(T)>::least;
};

template <std::integral T>
    requires(std::is_unsigned_v<T> && 2 * sizeof(T) <= sizeof(intmax_t))
struct double_integer<T>
{
    using type = boost::uint_t<CHAR_BIT * 2 * sizeof(T)>::least;
};

template <> struct double_integer<int64_t>
{
    using type = int128_t;
};

template <> struct double_integer<uint64_t>
{
    using type = uint128_t;
};

template <> struct double_integer<int128_t>
{
    using type = int256_t;
};

template <> struct double_integer<uint128_t>
{
    using type = uint256_t;
};

template <> struct double_integer<cpp_int>
{
    using type = cpp_int;
};

template <> struct double_integer<mpz_int>
{
    using type = mpz_int;
};

template <typename T> using double_integer_t = double_integer<T>::type;

template <std::integral T> struct half_integer
{
};

template <std::integral T>
    requires(std::is_signed_v<T>)
struct half_integer<T>
{
    using type = boost::int_t<CHAR_BIT / 2 * sizeof(T)>::least;
};

template <std::integral T>
    requires(std::is_unsigned_v<T>)
struct half_integer<T>
{
    using type = boost::uint_t<CHAR_BIT / 2 * sizeof(T)>::least;
};

template <typename T> using half_integer_t = half_integer<T>::type;

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
consteval int8_t operator""_i8(unsigned long long n) { return n; }
consteval int16_t operator""_i16(unsigned long long n) { return n; }
consteval int32_t operator""_i32(unsigned long long n) { return n; }
consteval int64_t operator""_i64(unsigned long long n) { return n; }
consteval int128_t operator""_i128(const char *str) { return fromString<int128_t>(str); }
consteval int256_t operator""_i256(const char *str) { return fromString<int256_t>(str); }
consteval int512_t operator""_i512(const char *str) { return fromString<int512_t>(str); }
consteval int1024_t operator""_i1024(const char *str) { return fromString<int1024_t>(str); }

consteval uint8_t operator""_u8(unsigned long long n) { return n; }
consteval uint16_t operator""_u16(unsigned long long n) { return n; }
consteval uint32_t operator""_u32(unsigned long long n) { return n; }
consteval uint64_t operator""_u64(unsigned long long n) { return n; }
consteval uint128_t operator""_u128(const char *str) { return fromString<uint128_t>(str); }
consteval uint256_t operator""_u256(const char *str) { return fromString<uint256_t>(str); }
consteval uint512_t operator""_u512(const char *str) { return fromString<uint512_t>(str); }
consteval uint1024_t operator""_u1024(const char *str) { return fromString<uint1024_t>(str); }

inline cpp_int operator""_cppi(const char *str) { return cpp_int{str}; }
inline mpz_int operator""_Z(const char *str) { return mpz_int{str}; }
inline mpq_rational operator""_Q(const char *str) { return mpq_rational{str}; }
inline mpf_float operator""_R(const char *str) { return mpf_float{str}; }
} // namespace literals

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
    if constexpr (requires(T a, T b) { a / b; })
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
        if (exponent & 1)
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
    using Td = double_integer_t<T>;
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
    if constexpr (!std::integral<T>)
        return sqrt(std::move(n));
    // boost::multiprecision::sqrt is constexpr, so take advantage of that in a constant-evaluated context.
    if (std::is_constant_evaluated())
        return boost::multiprecision::sqrt(std::move(n));
    // This constant is the first input where floor(sqrt(n)) returns the wrong value.
    if (n < 4'503'599'761'588'224)
        return sqrt(std::move(n));
    T x = sqrt(std::move(n));
    while (x * x > n)
        --x;
    return x;
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

/// GCD for int128.
inline int128_t gcd(int128_t m, int128_t n) noexcept { return boost::integer::gcd(m, n); }
/// LCM for int128.
inline int128_t lcm(int128_t m, int128_t n) noexcept { return boost::integer::lcm(m, n); }
/// GCD for uint128.
inline uint128_t gcd(uint128_t m, uint128_t n) noexcept { return boost::integer::gcd(m, n); }
/// LCM for uint128.
inline uint128_t lcm(uint128_t m, uint128_t n) noexcept { return boost::integer::lcm(m, n); }

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

/// A replacement for integer floor division.
template <typename T, typename U> constexpr auto fastDiv(T n, U d)
{
    auto res = T((double)n / (double)d);
    if (n < res * d)
    {
        --res;
        while (n < res * d)
            --res;
        return res;
    }
    if (n >= (res + 1) * d)
    {
        ++res;
        while (n >= (res + 1) * d)
            ++res;
        return res;
    }
    return res;
}

/// Wrapper for `mpz_divexact`.
inline mpz_int &divexact(mpz_int &dest, const mpz_int &a, const mpz_int &b)
{
    mpz_divexact(dest.backend().data(), a.backend().data(), b.backend().data());
    return dest;
}

/// Wrapper for `mpz_divexact_ui`.
inline mpz_int &divexact(mpz_int &dest, const mpz_int &a, uint32_t b)
{
    mpz_divexact_ui(dest.backend().data(), a.backend().data(), b);
    return dest;
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
    auto it = str.begin();

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

#ifdef BOOST_HAS_INT128
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
    constexpr auto operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) + std::forward<U>(b)))
        -> decltype(std::min(std::forward<T>(a), std::forward<U>(b)))
    {
        return std::min(std::forward<T>(a), std::forward<U>(b));
    }
};

struct maximum
{
    template <typename T, typename U>
    constexpr auto operator()(T &&a, U &&b) const noexcept(noexcept(std::forward<T>(a) + std::forward<U>(b)))
        -> decltype(std::max(std::forward<T>(a), std::forward<U>(b)))
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
    using reference = std::vector<T>::reference;
    using const_reference = std::vector<T>::const_reference;

    vector2d() = default;
    constexpr vector2d(size_t rows, size_t columns) : _rows(rows), _columns(columns), _data(rows * columns) {}
    constexpr vector2d(size_t rows, size_t columns, const T &value) : vector2d(rows, columns)
    {
        std::ranges::fill(_data, value);
    }

    /// Gets the linear index corresponding to a coordinate `(i, j)`.
    [[nodiscard]] constexpr int toIndex(int i, int j) const noexcept { return _columns * i + j; }

    /// `i`th row, `j`th column.
    constexpr const_reference operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr reference operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

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
            for (size_t j = 0; j < i; ++j)
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
            for (size_t j = 0; j < i; ++j)
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
} // namespace euler
