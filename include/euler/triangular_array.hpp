#pragma once

#include <array>
#include <cstddef>
#include <ostream>

namespace euler
{
#ifdef __cpp_multidimensional_subscript
// Triangular array. Consists of indices `(i, j)` where `j ≤ i`.
template <typename T, size_t N> class triangular_array
{
    std::array<T, N *(N + 1) / 2> _data;

    static constexpr size_t toIndex(size_t i, size_t j) noexcept
    {
        return j <= i ? (i * (i + 1) / 2) + j : (j * (j + 1) / 2) + i;
    }

  public:
    using iterator = std::array<T, N *(N + 1) / 2>::iterator;
    using const_iterator = std::array<T, N *(N + 1) / 2>::const_iterator;

    /// `i`th row, `j`th column.
    constexpr const T &operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr T &operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }
    [[nodiscard]] constexpr size_t height() const noexcept { return N; }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }

    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    [[nodiscard]] constexpr std::array<T, N *(N + 1) / 2> &data() noexcept { return _data; }
    [[nodiscard]] constexpr const std::array<T, N *(N + 1) / 2> &data() const noexcept { return _data; }

    /// `i`th row.
    constexpr std::span<T> operator[](size_t i) noexcept { return {_data.data() + i * (i + 1) / 2, i + 1}; }
    /// `i`th row.
    constexpr std::span<const T> operator[](size_t i) const noexcept { return {_data.data() + i * (i + 1) / 2, i + 1}; }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const triangular_array &v)
    {
        for (size_t i = 0; i < v.height(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                if (j != 0)
                    o << ", ";
                o << v[i, j];
            }
            o << '\n';
        }
        return o;
    }
};
#endif
} // namespace euler
