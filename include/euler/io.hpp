#pragma once

#include <chrono>
#include <fstream>

#include "decls.hpp"
#include "it/digits.hpp"

inline namespace euler
{
// namespace detail
// {
// /// Ignores ANSI colour sequences.
// inline size_t displayedWidth(std::string_view str)
// {
//     size_t result = 0;
//     const auto *p = str.begin();
//     for (; *p; ++p)
//     {
//         if (p[0] == '\e' && p[1] == '[')
//             while (*p != 'm')
//                 if (*p)
//                     ++p;
//                 else
//                     throw std::runtime_error("string terminates inside ANSI colour sequence");
//         else
//             ++result;
//     }
//     return result;
// }
// } // namespace detail

inline constexpr auto now = std::chrono::high_resolution_clock::now;

#ifdef BOOST_HAS_INT128
template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const int128_t &x)
{
    return o << boost::multiprecision::int128_t(x);
}

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const uint128_t &x)
{
    return o << boost::multiprecision::uint128_t(x);
}

template <typename CharT, typename Traits>
std::basic_istream<CharT, Traits> &operator>>(std::basic_istream<CharT, Traits> &is, int128_t &x)
{
    boost::multiprecision::int128_t a;
    is >> a;
    x = (int128_t)a;
    return is;
}

template <typename CharT, typename Traits>
std::basic_istream<CharT, Traits> &operator>>(std::basic_istream<CharT, Traits> &is, uint128_t &x)
{
    boost::multiprecision::int128_t a;
    is >> a;
    x = (uint128_t)a;
    return is;
}
#endif

template <typename CharT, typename Traits, typename T, typename U>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const std::pair<T, U> &v)
{
    std::basic_ostringstream<CharT, Traits> ss;
    ss.flags(o.flags());
    ss.imbue(o.getloc());
    ss.precision(o.precision());
    ss << '(' << v.first << ", " << v.second << ')';
    return o << std::move(ss).str();
}

template <typename CharT, typename Traits, typename... Args>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const std::tuple<Args...> &t)
{
    std::basic_ostringstream<CharT, Traits> ss;
    ss.flags(o.flags());
    ss.imbue(o.getloc());
    ss.precision(o.precision());
    bool first = true;
    ss << '(';
    apply([&](auto &&...args) { ((ss << (first ? "" : ", ") << args, first = false), ...); }, t);
    ss << ')';
    return o << std::move(ss).str();
}

template <typename CharT, typename Traits, typename T>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const std::optional<T> &x)
{
    if (x)
        return o << *x;
    return o << "none";
}

template <typename CharT, typename Traits, typename T>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, std::nullopt_t /*unused*/)
{
    return o << "none";
}

// Need to hide behind a namespace to avoid name collision with std::print.
namespace io
{
inline constexpr size_t defaultPrintLimit = 100;

/// Prints a non-range.
template <typename CharT = char, typename Traits = std::char_traits<CharT>, typename T>
std::basic_ostream<CharT, Traits> &
print(T &&x, size_t /* limit */ = defaultPrintLimit, std::basic_ostream<CharT, Traits> &o = std::cout,
      std::string_view /* delimiter */ = ", ", std::string_view /* open */ = "[", std::string_view /* close */ = "]")
{
    return o << std::forward<T>(x);
}

/// Prints a range.
template <typename CharT = char, typename Traits = std::char_traits<CharT>, std::ranges::range Range>
    requires(!is_string<Range>)
std::basic_ostream<CharT, Traits> &
print(Range &&r, size_t limit = defaultPrintLimit, std::basic_ostream<CharT, Traits> &o = std::cout,
      std::string_view delimiter = ", ", std::string_view open = "[", std::string_view close = "]")
{
    std::basic_ostringstream<CharT, Traits> ss;
    ss.flags(o.flags());
    ss.imbue(o.getloc());
    ss.precision(o.precision());
    ss << open;
    auto first = std::ranges::cbegin(r);
    auto last = std::ranges::cend(r);
    auto it = first;
    for (size_t i = 0; i < limit && it != last; ++i, ++it)
    {
        if (it != first)
            ss << delimiter;
        io::print(*it, limit, ss);
    }
    if (it != last)
    {
        if (it != first)
            ss << delimiter;
        ss << "...";
        if constexpr (std::ranges::sized_range<Range>)
            ss << " size=" << std::ranges::size(r);
    }
    ss << close;
    return o << std::move(ss).str();
}

/// Prints a duration.
template <typename CharT = char, typename Traits = std::char_traits<CharT>, typename Rep, typename Period>
std::basic_ostream<CharT, Traits> &
print(std::chrono::duration<Rep, Period> duration, size_t /*unused*/ = defaultPrintLimit,
      std::basic_ostream<CharT, Traits> &o = std::cout, std::string_view /*delimiter*/ = ", ",
      std::string_view /*open*/ = "[", std::string_view /*close*/ = "]")
{
    std::basic_ostringstream<CharT, Traits> ss;
    ss.flags(o.flags());
    ss.imbue(o.getloc());
    ss.precision(o.precision());
    auto dur = std::chrono::duration_cast<std::chrono::duration<double, std::nano>>(duration);
    auto count = dur.count();
    if (count < 1000)
        ss << count << " ns";
    else if (count < 1'000'000)
        ss << count / 1'000 << " us";
    else if (count < 1'000'000'000)
        ss << count / 1'000'000 << " ms";
    else
        ss << count / 1'000'000'000 << " s";
    return o << std::move(ss).str();
}

/// Prints anything and a newline.
template <typename CharT = char, typename Traits = std::char_traits<CharT>, typename T>
    requires(!is_string<T>)
std::basic_ostream<CharT, Traits> &
println(T &&x, size_t limit = defaultPrintLimit, std::basic_ostream<CharT, Traits> &o = std::cout,
        std::string_view delimiter = ", ", std::string_view open = "[", std::string_view close = "]")
{
    return print(std::forward<T>(x), limit, o, delimiter, open, close) << "\n";
}

/// Prints a range, each on its own line.
template <typename CharT = char, typename Traits = std::char_traits<CharT>, std::ranges::range Range>
    requires(!is_string<Range>)
std::basic_ostream<CharT, Traits> &printLines(Range &&r, size_t limit = defaultPrintLimit,
                                              std::basic_ostream<CharT, Traits> &o = std::cout)
{
    return println(std::forward<Range>(r), limit, o, "\n", "", "");
}
} // namespace io

template <typename CharT, typename Traits, std::ranges::range Range>
    requires(!is_string<Range>)
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Range &r)
{
    return io::print(r, io::defaultPrintLimit, o);
}

/// Converts any integral type to a string in the specified base.
template <integral2 T> constexpr std::string to_string(T n, int base = 10)
{
    bool neg = n < 0;
    if (neg)
        n = -n;
    static constexpr std::string_view digits =
        R"(0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz~!@#$%^&*()[]{}-+=;:'",./<>?\)";
    assert((size_t)base <= digits.size());
    std::string s(std::max(1UZ, countDigits(n, base) + neg), '0');
    if (neg)
        s[0] = '-';
    auto it = s.rbegin();
    it::digits(n, base)([&](int d) { *it++ = digits[d]; });
    return s;
}

/// Converts a duration to a friendly string.
template <typename Rep, typename Period> constexpr std::string to_string(const std::chrono::duration<Rep, Period> &d)
{
    std::ostringstream ss;
    io::print(d, io::defaultPrintLimit, ss);
    return std::move(ss).str();
}

/// Converts an integral number to its binary representation.
template <integral2 T> constexpr std::string bin(T n) { return to_string(std::move(n), 2); }

/// Reads the entire contents of a file as a string. Buggy with CRLF though.
inline std::string getFileContents(std::string_view fileName)
{
    if (std::ifstream fin(fileName.data(), std::ios::in); fin)
    {
        std::string result;
        fin.seekg(0, std::ios::end);
        result.resize(fin.tellg());
        fin.seekg(0, std::ios::beg);
        fin.read(result.data(), result.size());
        fin.close();
        return result;
    }
    throw(errno);
}

/// Measures the execution time of a given function and prints the timing information.
/// @param fn The function to be timed.
/// @param args The arguments to be passed to the function.
template <typename Callable, typename... Args>
    requires(!std::integral<Callable>)
void printTiming(Callable &&fn, Args &&...args)
{
    setConsoleToUtf8();
    std::ostringstream ss;
    auto const t1 = now();
    if constexpr (std::is_void_v<std::invoke_result_t<Callable, Args...>>)
    {
        std::forward<Callable>(fn)(std::forward<Args>(args)...);
        auto const elapsed = now() - t1;
        io::print(elapsed, io::defaultPrintLimit, ss);
        std::cout << std::move(ss).str() << '\n';
    }
    else
    {
        auto const result = std::forward<Callable>(fn)(std::forward<Args>(args)...);
        auto const elapsed = now() - t1;
        io::print(elapsed, io::defaultPrintLimit, ss);
        std::cout << std::move(ss).str() << " | " << result << '\n';
    }
}

/// Measures the execution time of a given function and prints the timing information.
/// @param repeat The number of times to repeat the function.
/// @param fn The function to be timed.
/// @param args The arguments to be passed to the function.
template <typename Callable, typename... Args> void printTiming(size_t repeat, Callable fn, Args... args)
{
    setConsoleToUtf8();
    std::ostringstream ss;
    auto t1 = now();
    if constexpr (std::is_void_v<std::invoke_result_t<Callable, Args...>>)
    {
        for (size_t i = 0; i < repeat; ++i)
            fn(args...);
        auto const elapsed = now() - t1;
        io::print(elapsed / (double)repeat, io::defaultPrintLimit, ss) << " (" << repeat << " iterations)";
        std::cout << std::move(ss).str() << '\n';
    }
    else
    {
        auto result = fn(args...);
        for (size_t i = 1; i < repeat; ++i)
            result = fn(args...);
        auto const elapsed = now() - t1;
        io::print(elapsed / (double)repeat, io::defaultPrintLimit, ss) << " (" << repeat << " iterations)";
        std::cout << std::move(ss).str() << " | " << result << '\n';
    }
}

/// Returns the string that would be output by `o << x`.
template <typename T, typename CharT, typename Traits>
std::string toStringWithFlags(const T &x, const std::basic_ostream<CharT, Traits> &o)
{
    std::basic_ostringstream<CharT, Traits> ss;
    ss.flags(o.flags());
    ss.imbue(o.getloc());
    ss.precision(o.precision());
    ss << x;
    return std::move(ss).str();
}

/// Prints a table. The row headers are `ir` and the column headers are `jr`.
/// @param ir The row headers.
/// @param jr The column headers.
/// @param f The function used to populate the cells of the table.
/// @param minCellWidth The minimum width of the table cells.
/// @param o The output stream.
template <std::ranges::range Range1, std::ranges::range Range2,
          std::invocable<std::ranges::range_value_t<Range1>, std::ranges::range_value_t<Range2>> Fun,
          typename CharT = char, typename Traits = std::char_traits<CharT>>
std::basic_ostream<CharT, Traits> &table(Range1 &&ir, Range2 &&jr, Fun f, int minCellWidth = 0,
                                         std::basic_ostream<CharT, Traits> &o = std::cout)
{
    int rowHeaderWidth = 0;
    std::vector<std::string> rowHeaders;
    for (auto &&i : ir)
    {
        auto &&s = rowHeaders.emplace_back(toStringWithFlags(i, o));
        rowHeaderWidth = std::max(rowHeaderWidth, (int)displayWidth(s));
    }
    std::vector<int> widths;
    for (auto &&j : jr)
        widths.push_back(std::max(minCellWidth, (int)displayWidth(toStringWithFlags(j, o))));
    std::vector<std::vector<std::string>> data;
    for (auto &&i : ir)
    {
        auto &v = data.emplace_back();
        size_t k = 0;
        for (auto &&j : jr)
        {
            const auto &s = v.emplace_back(toStringWithFlags(f(i, j), o));
            widths[k] = std::max(widths[k], (int)displayWidth(s));
            ++k;
        }
    }
    o << std::string(rowHeaderWidth + 1, ' ') << " │";
    auto wi = widths.begin();
    for (auto &&j : jr)
        o << std::setw(*wi++ + 1) << j;
    o << '\n';
    for (int i = 0; i <= rowHeaderWidth; ++i)
        o << "─";
    o << "─┼─";
    for (int w : widths)
        for (int i = 0; i <= w; ++i)
            o << "─";
    o << '\n';
    for (size_t i = 0; i < data.size(); ++i)
    {
        o << std::setw(rowHeaderWidth + rowHeaders[i].size() - displayWidth(rowHeaders[i]) + 1) << rowHeaders[i]
          << " │";
        for (size_t j = 0; j < data[i].size(); ++j)
            o << std::setw(widths[j] + data[i][j].size() - displayWidth(data[i][j]) + 1) << data[i][j];
        o << '\n';
    }
    return o;
}

/// Helpful using for legacy code.
using io::printLines;
using io::println;
} // namespace euler
