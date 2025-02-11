#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace ansi
{
/// Terminal color ansi escape codes.
enum style : uint8_t
{
    reset = 0,
    bold = 1,
    dim = 2,
    italic = 3,
    underline = 4,
    /// `slowBlink` and `rapidBlink` don't seem to work on modern terminals.
    slowBlink = 5,
    rapidBlink = 6,
    ///
    invert = 7,
    hide = 8,
    strike = 9,
    fontDefault = 10,
    doubleUnderline = 21,
    normalIntensity = 22,
    noItalic = 23,
    noUnderline = 24,
    noBlink = 25,
    /// Not known to be used in terminals.
    proportionalSpacing = 26,
    ///
    noReverse = 27,
    reveal = 28,
    noStrike = 29,
    black = 30,
    red = 31,
    green = 32,
    yellow = 33,
    blue = 34,
    magenta = 35,
    cyan = 36,
    white = 37,
    fgCustom = 38,
    fgDefault = 39,
    bgBlack = 40,
    bgRed = 41,
    bgGreen = 42,
    bgYellow = 43,
    bgBlue = 44,
    bgMagenta = 45,
    bgCyan = 46,
    bgWhite = 47,
    bgCustom = 48,
    bgDefault = 49,
    noProportionalSpacing = 50,
    frame = 51,
    encircle = 52,
    overline = 53,
    brightBlack = 90,
    brightRed = 91,
    brightGreen = 92,
    brightYellow = 93,
    brightBlue = 94,
    brightMagenta = 95,
    brightCyan = 96,
    brightWhite = 97,
    bgBrightBlack = 100,
    bgBrightRed = 101,
    bgBrightGreen = 102,
    bgBrightYellow = 103,
    bgBrightBlue = 104,
    bgBrightMagenta = 105,
    bgBrightCyan = 106,
    bgBrightWhite = 107
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
    for (uint8_t const x : rgb)
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
    for (uint8_t const x : rgb)
        s += ";" + std::to_string(x);
    s += "m";
    return s;
}

/// Returns the ansi code to set the background color to rgb.
inline std::string bg(uint8_t r, uint8_t g, uint8_t b) { return bg({r, g, b}); }

/// Clears the terminal screen.
constexpr std::string_view clear = "\033[2J";
} // namespace ansi
