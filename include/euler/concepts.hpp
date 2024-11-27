#pragma once

#include <boost/multiprecision/detail/number_base.hpp>
#include <execution>
#include <optional>

inline namespace euler
{
template <typename T>
concept integral2 = boost::multiprecision::number_category<T>::value ==
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
