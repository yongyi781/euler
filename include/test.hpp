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
    throw std::logic_error(std::string{message});
}

template <typename T, typename U, typename V = std::string_view>
void assertEqual(const T &actual, const U &expected, const V &name = {})
{
    using Tp = std::common_type_t<T, U>;
    using std::to_string;

    if (Tp(actual) != Tp(expected))
    {
        std::ostringstream ostr;
        if (name != V{})
            ostr << name << ": ";
        ostr << "Expected " << std::boolalpha;
        ostr << expected << ", got " << actual;
        fail(ostr.view());
    }
}
