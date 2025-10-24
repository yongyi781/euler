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
    throw std::logic_error(std::string{message});
}

template <typename T, typename U, typename V = std::string_view>
void assertEqual(const T &actual, const U &expected, const V &name = {})
{
    if (actual != expected)
    {
        std::ostringstream ostr;
        if (name != V{})
            ostr << name << ": ";
        ostr << "Expected " << std::boolalpha << expected << ", got " << actual;
        fail(ostr.view());
    }
}
