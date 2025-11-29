#pragma once

#include "ansi.hpp"
#include "euler/io.hpp"

#include <iostream>
#include <sstream>
#include <string>

inline void pass(std::string_view message)
{
    std::cerr << ansi::brightGreen << ansi::bold << "[PASS] " << ansi::reset << message << '\n';
}

[[noreturn]] inline void fail(std::string_view message)
{
    std::cerr << ansi::brightRed << ansi::bold << "[FAIL] " << ansi::reset << message << '\n';
    throw std::logic_error("Test failed");
}

template <typename T, typename U, typename V = std::string_view>
void assertEqual(const T &actual, const U &expected, const V &name = {},
                 const std::source_location &loc = std::source_location::current())
{
    if (actual != expected)
    {
        std::ostringstream ostr;
        if (!std::same_as<V, std::string_view> || name != V{})
            ostr << name << ": ";
        ostr << "Expected " << std::boolalpha << expected << ", got " << actual << " at " << loc.file_name() << ':'
             << loc.line();
        fail(ostr.view());
    }
}

template <typename T, typename U, typename V = std::string_view>
    requires std::floating_point<std::common_type_t<T, U>>
void assertNear(const T &actual, const U &expected, std::common_type_t<T, U> abs_tol = 1e-9,
                std::common_type_t<T, U> rel_tol = std::common_type_t<T, U>{0}, const V &name = {},
                const std::source_location &loc = std::source_location::current())
{
    using C = std::common_type_t<T, U>;

    const C a = static_cast<C>(actual);
    const C e = static_cast<C>(expected);

    // NaN handling: treat NaN as unequal unless both are NaN
    if (std::isnan(a) || std::isnan(e))
    {
        if (std::isnan(a) && std::isnan(e))
            return; // consider equal
        std::ostringstream ostr;
        if (!std::same_as<V, std::string_view> || name != V{})
            ostr << name << ": ";
        ostr << "Expected " << e << ", got " << a << " (NaN mismatch) at " << loc.file_name() << ':' << loc.line();
        fail(ostr.view());
    }

    // Infinity handling: must match exactly
    if (std::isinf(a) || std::isinf(e))
    {
        if (a == e)
            return;
        std::ostringstream ostr;
        if (!std::same_as<V, std::string_view> || name != V{})
            ostr << name << ": ";
        ostr << "Expected " << e << ", got " << a << " (infinity mismatch) at " << loc.file_name() << ':' << loc.line();
        fail(ostr.view());
    }

    const C diff = std::abs(a - e);
    const C scale = std::max(std::abs(a), std::abs(e));
    const C tol = std::max(abs_tol, rel_tol * scale);

    if (diff > tol)
    {
        std::ostringstream ostr;
        if (!std::same_as<V, std::string_view> || name != V{})
            ostr << name << ": ";
        ostr << "Expected " << e << " ± " << tol << " (abs_tol=" << abs_tol << ", rel_tol=" << rel_tol << "), got " << a
             << " (diff=" << diff << ") at " << loc.file_name() << ':' << loc.line();
        fail(ostr.view());
    }
}
