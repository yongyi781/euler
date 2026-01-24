#pragma once

#include "io.hpp"

#include <bit>
#include <cassert>
#include <ostream>
#include <vector>

namespace euler
{
/// Fenwick tree, also known as binary indexed tree.
template <typename T> class FenwickTree
{
    std::vector<T> data_;

    [[nodiscard]] static constexpr size_t parent(size_t i) noexcept { return i - (i & -i); }
    [[nodiscard]] static constexpr size_t next(size_t i) noexcept { return i + (i & -i); }

  public:
    /// Creates an empty Fenwick tree.
    FenwickTree() = default;
    /// Constructs a new Fenwick tree with the given size.
    constexpr explicit FenwickTree(size_t n) : data_(n + 1, T{}) {}

    /// Constructs a new Fenwick tree with the given size and initial value.
    constexpr FenwickTree(size_t n, T value) : data_(n + 1)
    {
        for (size_t i = 1; i <= n; i++)
        {
            data_[i] += value;
            if (next(i) <= n)
                data_[next(i)] += data_[i];
        }
    }

    /// Constructs a new Fenwick tree with the given size populated by the given function.
    template <std::invocable<size_t> Fun> constexpr FenwickTree(size_t n, Fun fn) : data_(n + 1)
    {
        for (size_t i = 1; i <= n; i++)
        {
            data_[i] += fn(i - 1);
            if (next(i) <= n)
                data_[next(i)] += data_[i];
        }
    }

    /// Constructs a new Fenwick tree with the given initial data.
    template <std::forward_iterator It> constexpr FenwickTree(It first, It last) : data_(std::distance(first, last) + 1)
    {
        size_t const n = size();
        for (size_t i = 1; i <= n; i++)
        {
            data_[i] += *first++;
            if (next(i) <= n)
                data_[next(i)] += data_[i];
        }
    }

    std::vector<T> &data() noexcept { return data_; }
    [[nodiscard]] constexpr const std::vector<T> &data() const noexcept { return data_; }

    /// Returns the logical size of the Fenwick tree. It is 1 less than the physical size.
    [[nodiscard]] constexpr size_t size() const noexcept { return data_.size() - 1; }

    /// Resizes the Fenwick tree.
    constexpr void resize(size_t newSize) { data_.resize(newSize + 1); }

    /// Returns the sum of the elements up to the given index in the Fenwick tree. O(log n).
    [[nodiscard]] constexpr T sum(size_t i) const
    {
        assert(i < size());
        T res{};
        for (++i; i != 0; i = parent(i))
            res += data_[i];
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
        for (++i; i < data_.size(); i = next(i))
            data_[i] += v;
    }

    /// Returns the value at a specific index. O(log n).
    constexpr T operator[](size_t i) const { return sum(i, i); }

    // Returns the first index `i` such that `sum(i) ≥ k`. Requires all the values to be non-negative for monotonicity.
    [[nodiscard]] constexpr size_t first_ge(T k) const
    {
        size_t i = 0;
        T s{};
        for (size_t mask = std::bit_floor(size()); mask > 0; mask >>= 1)
        {
            if (i + mask <= size() && s + data_[i + mask] < k)
            {
                i += mask;
                s += data_[i];
            }
        }
        return i;
    }

    // Returns the first index `i` such that `sum(i) > k`. Requires all the values to be non-negative for monotonicity.
    [[nodiscard]] constexpr size_t first_gt(T k) const
    {
        size_t i = 0;
        T s = 0;
        for (size_t mask = std::bit_floor(size()); mask > 0; mask >>= 1)
        {
            if (i + mask <= size() && s + data_[i + mask] <= k)
            {
                i += mask;
                s += data_[i];
            }
        }
        return i;
    }

    /// Converts this Fenwick tree to a `std::vector` of prefix sums.
    [[nodiscard]] std::vector<T> toPrefixSum() const
    {
        size_t const n = size();
        std::vector<T> res(n);
        for (size_t i = 1; i <= n; ++i)
        {
            size_t const j = parent(i);
            res[i - 1] = data_[i] + (j == 0 ? T{} : res[j - 1]);
        }
        return res;
    }

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
