#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace ansi
{
/// Terminal color ansi escape codes.
enum style : uint8_t
{
    reset,
    bold,
    dim,
    italic,
    underline,
    /// `slowBlink` and `rapidBlink` don't seem to work on modern terminals.
    slowBlink,
    rapidBlink,
    ///
    invert,
    hide,
    strike,
    fontDefault,
    doubleUnderline = 21,
    normalIntensity,
    noItalic,
    noUnderline,
    noBlink,
    /// Not known to be used in terminals.
    proportionalSpacing,
    ///
    noReverse,
    reveal,
    noStrike,
    black,
    red,
    green,
    yellow,
    blue,
    magenta,
    cyan,
    white,
    fgCustom,
    fgDefault,
    bgBlack,
    bgRed,
    bgGreen,
    bgYellow,
    bgBlue,
    bgMagenta,
    bgCyan,
    bgWhite,
    bgCustom,
    bgDefault,
    noProportionalSpacing,
    frame,
    encircle,
    overline,
    brightBlack = 90,
    brightRed,
    brightGreen,
    brightYellow,
    brightBlue,
    brightMagenta,
    brightCyan,
    brightWhite,
    bgBrightBlack = 100,
    bgBrightRed,
    bgBrightGreen,
    bgBrightYellow,
    bgBrightBlue,
    bgBrightMagenta,
    bgBrightCyan,
    bgBrightWhite
};

/// Converts an ANSI style to a string.
inline std::string str(style n) { return "\033[" + std::to_string((int)n) + "m"; }

/// Prints an ansi style code to a stream.
template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &os, style n)
{
    return os << str(n);
}

/// Returns the ansi code to set the foreground color to rgb.
inline std::string fg(std::array<uint8_t, 3> rgb)
{
    std::string s = "\033[38;2";
    for (uint8_t x : rgb)
        s += ";" + std::to_string(x);
    s += "m";
    return s;
}

/// Returns the ansi code to set the foreground color to rgb.
inline std::string fg(uint8_t r, uint8_t g, uint8_t b) { return fg({r, g, b}); }

/// Returns the ansi code to set the background color to rgb.
inline std::string bg(std::array<uint8_t, 3> rgb)
{
    std::string s = "\033[48;2";
    for (uint8_t x : rgb)
        s += ";" + std::to_string(x);
    s += "m";
    return s;
}

/// Returns the ansi code to set the background color to rgb.
inline std::string bg(uint8_t r, uint8_t g, uint8_t b) { return bg({r, g, b}); }

/// Clears the terminal screen.
constexpr std::string_view clear = "\033[2J";
} // namespace ansi
