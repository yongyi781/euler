#pragma once

#include <cstddef>
#include <ostream>
#include <vector>

inline namespace euler
{
#ifdef __cpp_multidimensional_subscript
// Triangular vector. Consists of indices `(i, j)` where `j ≤ i`.
template <typename T> class triangular_vector
{
  public:
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

    triangular_vector() = default;
    constexpr triangular_vector(size_t height) : _height(height), _data(toIndex(height, 0)) {}

    /// `i`th row, `j`th column.
    constexpr const T &operator[](size_t i, size_t j) const noexcept { return _data[toIndex(i, j)]; }
    /// `i`th row, `j`th column.
    constexpr T &operator[](size_t i, size_t j) noexcept { return _data[toIndex(i, j)]; }

    [[nodiscard]] constexpr size_t height() const noexcept { return _height; }

    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }

    [[nodiscard]] constexpr iterator begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return _data.begin(); }

    [[nodiscard]] constexpr iterator end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return _data.end(); }

    constexpr void resize(size_t newHeight)
    {
        _height = newHeight;
        _data.resize(toIndex(newHeight, 0));
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const triangular_vector &v)
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

  private:
    size_t _height = 0;
    std::vector<T> _data;

    static constexpr size_t toIndex(size_t i, size_t j) noexcept { return (i * (i + 1) / 2) + j; }
};
#endif
} // namespace euler
