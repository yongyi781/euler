#pragma once

#include <cassert>
#include <ostream>
#include <vector>

inline namespace euler
{
template <typename T> class FenwickTree
{
  public:
    /// Constructs a new Fenwick tree with the given size.
    explicit constexpr FenwickTree(size_t n) : _data(n + 1, T{}) {}

    /// Constructs a new Fenwick tree with the given size and initial value.
    constexpr FenwickTree(size_t n, T value) : _data(n + 1)
    {
        for (size_t i = 1; i <= n; i++)
        {
            _data[i] += value;
            if (getNext(i) <= n)
                _data[getNext(i)] += _data[i];
        }
    }

    /// Constructs a new Fenwick tree with the given size populated by the given function.
    template <std::invocable<size_t> Fun> constexpr FenwickTree(size_t n, Fun fn) : _data(n + 1)
    {
        for (size_t i = 1; i <= n; i++)
        {
            _data[i] += fn(i - 1);
            if (getNext(i) <= n)
                _data[getNext(i)] += _data[i];
        }
    }

    /// Constructs a new Fenwick tree with the given initial data.
    template <std::forward_iterator It> constexpr FenwickTree(It first, It last) : _data(std::distance(first, last) + 1)
    {
        size_t const n = size();
        for (size_t i = 1; i <= n; i++)
        {
            _data[i] += *first++;
            if (getNext(i) <= n)
                _data[getNext(i)] += _data[i];
        }
    }

    std::vector<T> &data() { return _data; }
    [[nodiscard]] constexpr const std::vector<T> &data() const { return _data; }

    /// Returns the size of the data.
    [[nodiscard]] constexpr size_t size() const { return _data.size() - 1; }

    /// Returns the sum of the elements up to the given index in the Fenwick tree.
    ///
    /// @param i the index up to which the sum should be calculated.
    [[nodiscard]] constexpr T getSum(size_t i) const
    {
        assert(i < size());
        ++i;
        T res{};
        while (i > 0)
        {
            res += _data[i];
            i = getParent(i);
        }
        return res;
    }

    /// Calculates the sum of elements in a range.
    ///
    /// @param i The starting index of the range.
    /// @param j The ending index of the range.
    ///
    /// @return The sum of elements in the range [i, j].
    [[nodiscard]] constexpr T getSum(size_t i, size_t j) const
    {
        assert(i <= j && j < size());
        return getSum(j) - (i == 0 ? T{} : getSum(i - 1));
    }

    /// Updates the Fenwick tree by adding the value `v` to position `i`.
    ///
    /// @param i the starting index from which to update the elements
    /// @param v the value to add.
    constexpr void update(size_t i, T v)
    {
        ++i;
        while (i < _data.size())
        {
            _data[i] += v;
            i = getNext(i);
        }
    }

    /// Indexing operator to access the value at a specific index. O(log n) time.
    ///
    /// @param i the index to access.
    /// @return the value at the specified index.
    constexpr T operator[](size_t i) const { return getSum(i, i); }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const FenwickTree<T> &ft)
    {
        return o << ft.data();
    }

  private:
    std::vector<T> _data;

    [[nodiscard]] static constexpr size_t getParent(size_t i) { return i - (i & (-i)); }
    [[nodiscard]] static constexpr size_t getNext(size_t i) { return i + (i & (-i)); }
};
} // namespace euler
