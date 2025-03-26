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
inline void fail(std::string_view message)
{
    std::cerr << ansi::brightRed << ansi::bold << "[FAIL] " << ansi::reset << message << '\n';
}

template <typename T, typename U> void assertEqual(const T &actual, const U &expected)
{
    using Tp = std::common_type_t<T, U>;
    using std::to_string;

    if (Tp(actual) != Tp(expected))
    {
        std::ostringstream ostr;
        ostr << "Expected ";
        if constexpr (std::same_as<Tp, bool>)
            ostr << (expected ? "true, got false" : "false, got true");
        else
            ostr << expected << ", got " << actual;
        fail(ostr.view());
        throw std::logic_error("Test failed");
    }
}
