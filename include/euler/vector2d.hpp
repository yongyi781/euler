#pragma once

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>

namespace euler
{
#ifdef __cpp_multidimensional_subscript
/// 2D vector.
template <typename T> class vector2d
{
    size_t _dim0 = 0;
    size_t _dim1 = 0;
    std::vector<T> _data;

    /// Gets the linear index corresponding to a coordinate `(i, j)`.
    [[nodiscard]] constexpr size_t toIndex(size_t i, size_t j) const noexcept { return _dim1 * i + j; }

  public:
    using value_type = T;
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;
    using reference = std::vector<T>::reference;
    using const_reference = std::vector<T>::const_reference;

    vector2d() = default;
    constexpr vector2d(size_t dim0, size_t dim1) : _dim0(dim0), _dim1(dim1), _data(dim0 * dim1) {}
    constexpr vector2d(size_t dim0, size_t dim1, const T &value) : vector2d(dim0, dim1)
    {
        std::ranges::fill(_data, value);
    }

    [[nodiscard]] constexpr std::array<size_t, 2> dims() const noexcept { return {_dim0, _dim1}; }
    [[nodiscard]] constexpr size_t dim(size_t d) const noexcept { return dims()[d]; }
    [[nodiscard]] constexpr size_t dim0() const noexcept { return _dim0; }
    [[nodiscard]] constexpr size_t dim1() const noexcept { return _dim1; }

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
    constexpr std::span<T> operator[](size_t i) noexcept { return {_data.data() + i * _dim1, _dim1}; }
    /// `i`th row.
    constexpr std::span<const T> operator[](size_t i) const noexcept { return {_data.data() + i * _dim1, _dim1}; }

    constexpr auto operator<=>(const vector2d &other) const = default;

    /// Zeros out the vector. This method does not set size to 0.
    constexpr void clear() noexcept { std::ranges::fill(_data, T{}); }

    /// Resizes (and potentially reshapes) the vector.
    constexpr void resize(size_t dim0, size_t dim1)
    {
        _dim0 = dim0;
        _dim1 = dim1;
        _data.resize(dim1 * dim0);
    }

    constexpr vector2d transpose() const
    {
        vector2d res(_dim1, _dim0);
        for (size_t i = 0; i < _dim0; ++i)
            for (size_t j = 0; j < _dim1; ++j)
                res[j, i] = (*this)[i, j];
        return res;
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
        ss << v._dim0 << "×" << v._dim1 << " vector2d:\n";
        for (size_t i = 0; i < v.dim0(); ++i)
        {
            for (size_t j = 0; j < v.dim1(); ++j)
                ss << std::setw(maxWidth + 1) << v[i, j];
            ss << '\n';
        }
        return o << std::move(ss).str();
    }
};
#endif
} // namespace euler
