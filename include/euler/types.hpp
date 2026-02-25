#pragma once

#include <boost/multiprecision/integer.hpp>
#include <execution>
#include <optional>

// Helpers for std::is_signed and std::make_signed with boost integer types.
template <std::size_t MinBits, std::size_t MaxBits, boost::multiprecision::cpp_integer_type SignType,
          boost::multiprecision::cpp_int_check_type Checked, class Allocator>
struct std::make_signed<boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>>
{
    using type = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
        MinBits, MaxBits, boost::multiprecision::signed_magnitude, Checked, Allocator>>;
};

template <std::size_t MinBits, std::size_t MaxBits, boost::multiprecision::cpp_integer_type SignType,
          boost::multiprecision::cpp_int_check_type Checked, class Allocator>
struct std::is_signed<boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>>
    : public integral_constant<bool, SignType == boost::multiprecision::signed_magnitude ||
                                         SignType == boost::multiprecision::signed_packed>
{
};

// Other backends are signed as far as I can tell.
template <typename Backend, boost::multiprecision::expression_template_option ExpressionTemplates>
struct std::make_signed<boost::multiprecision::number<Backend, ExpressionTemplates>>
{
    using type = boost::multiprecision::number<Backend, ExpressionTemplates>;
};

template <typename Backend, boost::multiprecision::expression_template_option ExpressionTemplates>
struct std::is_signed<boost::multiprecision::number<Backend, ExpressionTemplates>> : public true_type
{
};

namespace euler
{
#ifdef __SIZEOF_INT128__
using int128_t = __int128;
using uint128_t = unsigned __int128;
#else
using int128_t = boost::multiprecision::int128_t;
using uint128_t = boost::multiprecision::uint128_t;
#endif
using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;

using boost::multiprecision::mpf_float;
using boost::multiprecision::mpq_rational;
using boost::multiprecision::mpz_int;

// Rust-style type abbreviations
using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;
using i128 = int128_t;

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using u128 = uint128_t;

template <typename T> using PrimePower = std::pair<T, int>;

template <typename T> struct double_integer
{
};

template <std::signed_integral T>
    requires(2 * sizeof(T) <= sizeof(intmax_t))
struct double_integer<T>
{
    using type = boost::int_t<CHAR_BIT * 2 * sizeof(T)>::least;
};

template <std::unsigned_integral T>
    requires(2 * sizeof(T) <= sizeof(intmax_t))
struct double_integer<T>
{
    using type = boost::uint_t<CHAR_BIT * 2 * sizeof(T)>::least;
};

template <> struct double_integer<i64>
{
    using type = i128;
};

template <> struct double_integer<u64>
{
    using type = u128;
};

template <std::size_t MinBits, std::size_t MaxBits, boost::multiprecision::cpp_integer_type SignType,
          boost::multiprecision::cpp_int_check_type Checked, class Allocator>
struct double_integer<boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>>
{
    using type = boost::multiprecision::number<
        boost::multiprecision::cpp_int_backend<2 * MinBits, 2 * MaxBits, SignType, Checked, Allocator>>;
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

template <typename T> struct half_integer
{
    using type = T;
};

template <std::signed_integral T> struct half_integer<T>
{
    using type = boost::int_t<CHAR_BIT / 2 * sizeof(T)>::least;
};

template <std::unsigned_integral T> struct half_integer<T>
{
    using type = boost::uint_t<CHAR_BIT / 2 * sizeof(T)>::least;
};

template <typename T> using half_integer_t = half_integer<T>::type;

template <typename T, size_t D0, size_t D1> using array2d = std::array<std::array<T, D1>, D0>;

template <typename T, size_t D0, size_t D1, size_t D2>
using array3d = std::array<std::array<std::array<T, D2>, D1>, D0>;

// ==== Concepts ==============================================================

template <typename T>
concept integral2 = std::numeric_limits<T>::is_integer;

template <typename T>
concept is_string =
    std::same_as<std::decay_t<T>, const char *> || std::same_as<std::decay_t<T>, char *> ||
    std::same_as<std::decay_t<T>, const wchar_t *> || std::same_as<std::decay_t<T>, wchar_t *> ||
    std::same_as<std::decay_t<T>, std::string> || std::same_as<std::decay_t<T>, std::wstring> ||
    std::same_as<std::decay_t<T>, std::u16string> || std::same_as<std::decay_t<T>, std::u32string> ||
    std::same_as<std::decay_t<T>, std::u8string> || std::same_as<std::decay_t<T>, std::string_view> ||
    std::same_as<std::decay_t<T>, std::wstring_view> || std::same_as<std::decay_t<T>, std::u16string_view> ||
    std::same_as<std::decay_t<T>, std::u32string_view> || std::same_as<std::decay_t<T>, std::u8string_view>;

template <typename T>
concept is_optional = std::same_as<T, std::optional<typename T::value_type>>;

template <typename T>
concept execution_policy = std::is_execution_policy_v<std::decay_t<T>>;

template <typename T>
concept parallel_policy =
    execution_policy<T> && (std::is_same_v<std::decay_t<T>, std::execution::parallel_policy> ||
                            std::is_same_v<std::decay_t<T>, std::execution::parallel_unsequenced_policy>);
} // namespace euler
