#pragma once

#include <cassert>
#include <ostream>
#include <vector>

#include "io.hpp"

inline namespace euler
{
/// Fenwick tree, also known as binary indexed tree.
template <typename T> class FenwickTree
{
    std::vector<T> _data;

    [[nodiscard]] static constexpr size_t parent(size_t i) noexcept { return i - (i & -i); }
    [[nodiscard]] static constexpr size_t next(size_t i) noexcept { return i + (i & -i); }

  public:
    /// Creates an empty Fenwick tree.
    FenwickTree() = default;
    /// Constructs a new Fenwick tree with the given size.
    constexpr explicit FenwickTree(size_t n) : _data(n + 1, T{}) {}

    /// Constructs a new Fenwick tree with the given size and initial value.
    constexpr FenwickTree(size_t n, T value) : _data(n + 1)
    {
        for (size_t i = 1; i <= n; i++)
        {
            _data[i] += value;
            if (next(i) <= n)
                _data[next(i)] += _data[i];
        }
    }

    /// Constructs a new Fenwick tree with the given size populated by the given function.
    template <std::invocable<size_t> Fun> constexpr FenwickTree(size_t n, Fun fn) : _data(n + 1)
    {
        for (size_t i = 1; i <= n; i++)
        {
            _data[i] += fn(i - 1);
            if (next(i) <= n)
                _data[next(i)] += _data[i];
        }
    }

    /// Constructs a new Fenwick tree with the given initial data.
    template <std::forward_iterator It> constexpr FenwickTree(It first, It last) : _data(std::distance(first, last) + 1)
    {
        size_t const n = size();
        for (size_t i = 1; i <= n; i++)
        {
            _data[i] += *first++;
            if (next(i) <= n)
                _data[next(i)] += _data[i];
        }
    }

    std::vector<T> &data() noexcept { return _data; }
    [[nodiscard]] constexpr const std::vector<T> &data() const noexcept { return _data; }

    /// Returns the size of the Fenwick tree.
    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size() - 1; }

    /// Resizes the Fenwick tree.
    constexpr void resize(size_t newSize) { _data.resize(newSize + 1); }

    /// Returns the sum of the elements up to the given index in the Fenwick tree. O(log n).
    [[nodiscard]] constexpr T sum(size_t i) const
    {
        assert(i < size());
        ++i;
        T res{};
        while (i != 0)
        {
            res += _data[i];
            i = parent(i);
        }
        return res;
    }

    /// Calculates the sum of elements in the range [i, j]. O(2 * log n).
    [[nodiscard]] constexpr T sum(size_t i, size_t j) const
    {
        assert(j < size());
        if (i > j)
            return T{};
        return sum(j) - (i == 0 ? T{} : sum(i - 1));
    }

    /// Adds the value `v` to position `i`. O(log n).
    constexpr void add(size_t i, T v)
    {
        ++i;
        while (i < _data.size())
        {
            _data[i] += v;
            i = next(i);
        }
    }

    /// Returns the value at a specific index. O(log n).
    constexpr T operator[](size_t i) const { return sum(i, i); }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const FenwickTree<T> &t)
    {
        o << "Fenwick tree:\n  Values: [";
        for (size_t i = 0; i < t.size(); i++)
            o << (i == 0 ? "" : ", ") << t[i];
        o << "]\n  Raw data: " << t.data() << '\n';
        return o;
    }
};
} // namespace euler
