#pragma once

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>

inline namespace euler
{
#ifdef __cpp_multidimensional_subscript
/// 2D vector.
template <typename T> class vector2d
{
  public:
    using value_type = T;
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;
    using reference = std::vector<T>::reference;
    using const_reference = std::vector<T>::const_reference;

    vector2d() = default;
    constexpr vector2d(size_t rows, size_t columns) : _rows(rows), _columns(columns), _data(rows * columns) {}
    constexpr vector2d(size_t rows, size_t columns, const T &value) : vector2d(rows, columns)
    {
        std::ranges::fill(_data, value);
    }

    [[nodiscard]] constexpr size_t rows() const noexcept { return _rows; }
    [[nodiscard]] constexpr size_t columns() const noexcept { return _columns; }

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
    constexpr std::span<T> operator[](size_t i) noexcept { return {_data.data() + i * _columns, _columns}; }
    /// `i`th row.
    constexpr std::span<const T> operator[](size_t i) const noexcept { return {_data.data() + i * _columns, _columns}; }

    /// Zeros out the vector. This method does not set size to 0.
    constexpr void clear() noexcept { std::ranges::fill(_data, T{}); }

    constexpr void resize(size_t rows, size_t columns)
    {
        _rows = rows;
        _columns = columns;
        _data.resize(columns * rows);
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const vector2d &v)
    {
        int maxWidth = o.width();
        for (auto &&x : v._data)
            maxWidth = std::max(maxWidth, (int)toStringWithFlags(x, o).size());
        std::basic_ostringstream<CharT, Traits> ss;
        ss.flags(o.flags());
        ss.imbue(o.getloc());
        ss.precision(o.precision());
        ss << v._columns << "×" << v._rows << " vector2d:\n";
        for (size_t i = 0; i < v.rows(); ++i)
        {
            for (size_t j = 0; j < v.columns(); ++j)
                ss << std::setw(maxWidth + 1) << v[i, j];
            ss << '\n';
        }
        return o << std::move(ss).str();
    }

  private:
    size_t _rows = 0;
    size_t _columns = 0;
    std::vector<T> _data;

    /// Gets the linear index corresponding to a coordinate `(i, j)`.
    [[nodiscard]] constexpr size_t toIndex(size_t i, size_t j) const noexcept { return _columns * i + j; }
};
#endif
} // namespace euler
