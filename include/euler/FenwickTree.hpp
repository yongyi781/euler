#pragma once

#include <cassert>
#include <vector>

inline namespace euler
{
template <typename T> class FenwickTree
{
  public:
    explicit FenwickTree(std::size_t n) : _data(n + 1, T()) {}

    std::vector<T> &data() { return _data; }

    /**
     * Returns the size of the data.
     *
     * @return the size of the data.
     */
    [[nodiscard]] std::size_t size() const { return _data.size(); }

    /**
     * Returns the sum of the elements up to the given index in the Fenwick tree.
     *
     * @param i the index up to which the sum should be calculated.
     *
     * @return the sum of the elements up to the given index.
     */
    [[nodiscard]] T getSum(std::size_t i) const
    {
        assert(i < _data.size());
        T sum = 0;
        ++i;
        while (i > 0)
        {
            sum += _data[i];
            i = getParent(i);
        }
        return sum;
    }

    /**
     * Calculates the sum of elements in a range.
     *
     * @param i The starting index of the range.
     * @param j The ending index of the range.
     *
     * @return The sum of elements in the range [i, j].
     */
    [[nodiscard]] T getSum(size_t i, size_t j) const
    {
        assert(i <= j);
        return getSum(j) - (i == 0 ? T() : getSum(i - 1));
    }

    /**
     * Updates the Fenwick tree by adding the value `v` to position `i`.
     *
     * @param i the starting index from which to update the elements
     * @param v the value to add.
     */
    void update(size_t i, T v)
    {
        ++i;
        while (i < _data.size())
        {
            _data[i] += v;
            i = getNext(i);
        }
    }

  private:
    std::vector<T> _data;

    [[nodiscard]] size_t getParent(size_t i) const { return i - (i & (-i)); }
    [[nodiscard]] size_t getNext(size_t i) const { return i + (i & (-i)); }
};
} // namespace euler
