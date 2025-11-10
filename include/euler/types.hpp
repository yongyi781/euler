#pragma once

#include <boost/multiprecision/integer.hpp>
#include <execution>
#include <optional>

namespace std
{
// Helpers for std::is_signed and std::make_signed with boost integer types.
template <std::size_t MinBits, std::size_t MaxBits, boost::multiprecision::cpp_integer_type SignType,
          boost::multiprecision::cpp_int_check_type Checked, class Allocator>
struct make_signed<boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>>
{
    using type = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
        MinBits, MaxBits, boost::multiprecision::signed_magnitude, Checked, Allocator>>;
};

template <std::size_t MinBits, std::size_t MaxBits, boost::multiprecision::cpp_integer_type SignType,
          boost::multiprecision::cpp_int_check_type Checked, class Allocator>
struct is_signed<boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>>
    : public integral_constant<bool, SignType == boost::multiprecision::signed_magnitude ||
                                         SignType == boost::multiprecision::signed_packed>
{
};

// Other backends are signed as far as I can tell.
template <typename Backend, boost::multiprecision::expression_template_option ExpressionTemplates>
struct make_signed<boost::multiprecision::number<Backend, ExpressionTemplates>>
{
    using type = boost::multiprecision::number<Backend, ExpressionTemplates>;
};

template <typename Backend, boost::multiprecision::expression_template_option ExpressionTemplates>
struct is_signed<boost::multiprecision::number<Backend, ExpressionTemplates>> : public true_type
{
};
} // namespace std

#ifdef __SIZEOF_INT128__
using int128_t = __int128;
using uint128_t = unsigned __int128;
#else
using int128_t = boost::multiprecision::int128_t;
using uint128_t = boost::multiprecision::uint128_t;
#endif
using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using boost::multiprecision::int256_t;
using boost::multiprecision::int512_t;

using boost::multiprecision::int1024_t;
using int2048_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    2048, 2048, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using int4096_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    4096, 4096, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using int8192_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    8192, 8192, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void>>;
using boost::container::static_vector;

using boost::multiprecision::mpf_float;
using boost::multiprecision::mpq_rational;
using boost::multiprecision::mpz_int;
using boost::multiprecision::uint1024_t;
using boost::multiprecision::uint256_t;
using boost::multiprecision::uint512_t;

namespace euler
{
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

template <typename T> struct half_integer
{
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

// ==== Concepts ==============================================================

template <typename T>
concept integral2 = boost::multiprecision::number_category<std::decay_t<T>>::value ==
                    boost::multiprecision::number_category_type::number_kind_integer;

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
} // namespace euler

// Rust-style type abbreviations
using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;
using i128 = int128_t;
using i256 = int256_t;
using i512 = int512_t;
using i1024 = int1024_t;

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using u128 = uint128_t;
using u256 = uint256_t;
using u512 = uint512_t;
using u1024 = uint1024_t;
