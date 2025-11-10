#pragma once

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <vector>

namespace euler
{
#ifdef __cpp_multidimensional_subscript
// Triangular vector. Consists of indices `(i, j)` where `j ≤ i`.
template <typename T> class triangular_vector
{
    size_t _height = 0;
    std::vector<T> _data;

    /// Gets the linear index corresponding to a coordinate `(i, j)`.
    [[nodiscard]] static constexpr size_t toIndex(size_t i, size_t j) noexcept { return (i * (i + 1) / 2) + j; }

  public:
    using value_type = T;
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;
    using reference = std::vector<T>::reference;
    using const_reference = std::vector<T>::const_reference;

    triangular_vector() = default;
    constexpr triangular_vector(size_t height) : _height(height), _data(toIndex(height, 0)) {}
    constexpr triangular_vector(size_t height, const T &value) : _height(height), _data(toIndex(height, 0))
    {
        std::ranges::fill(_data, value);
    }

    [[nodiscard]] constexpr size_t height() const noexcept { return _height; }

    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }

    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    [[nodiscard]] constexpr std::vector<T> &data() noexcept { return _data; }
    [[nodiscard]] constexpr const std::vector<T> &data() const noexcept { return _data; }

    /// `i`th row, `j`th column.
    constexpr const_reference operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr reference operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

    /// `i`th row.
    constexpr std::span<T> operator[](size_t i) noexcept { return {_data.data() + i * (i + 1) / 2, i + 1}; }
    /// `i`th row.
    constexpr std::span<const T> operator[](size_t i) const noexcept { return {_data.data() + i * (i + 1) / 2, i + 1}; }

    constexpr auto operator<=>(const triangular_vector &other) const = default;

    /// Zeros out the vector. This method does not set size to 0.
    constexpr void clear() noexcept { std::ranges::fill(_data, T{}); }

    constexpr void resize(size_t newHeight)
    {
        _height = newHeight;
        _data.resize(toIndex(newHeight, 0));
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const triangular_vector &v)
    {
        int maxWidth = o.width();
        for (auto &&x : v._data)
            maxWidth = std::max(maxWidth, (int)toStringWithFlags(x, o).size());
        std::basic_ostringstream<CharT, Traits> ss;
        ss.flags(o.flags());
        ss.imbue(o.getloc());
        ss.precision(o.precision());
        ss << "Triangular vector of height " << v.height() << ":\n";
        for (size_t i = 0; i < v.height(); ++i)
        {
            for (size_t j = 0; j <= i; ++j)
                ss << std::setw(maxWidth + 1) << v[i, j];
            ss << '\n';
        }
        return o << std::move(ss).str();
    }
};
#endif
} // namespace euler
