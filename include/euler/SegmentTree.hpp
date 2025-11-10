#pragma once

#include <ostream>
#include <ranges>
#include <vector>

namespace euler
{
/// A segment tree for an arbitrary monoid.
template <typename Node> class SegmentTree
{
    std::vector<Node> _data;

  public:
    /// Creates an empty segment tree of given size.
    constexpr SegmentTree(size_t n) : _data(2 * n) {}

    /// Creates an empty segment tree of given size.
    constexpr SegmentTree(size_t n, const Node &value) : _data(2 * n)
    {
        std::ranges::fill(_data.begin() + size(), _data.end(), value);
        for (size_t i = size() - 1; i > 0; --i)
            _data[i] = _data[2 * i] + _data[2 * i + 1];
    }

    /// Creates a segment tree from a range of elements.
    template <std::ranges::sized_range Range> constexpr SegmentTree(Range &&init) : SegmentTree(init.size())
    {
        std::ranges::move(std::forward<Range>(init), _data.begin() + size());
        for (size_t i = size() - 1; i > 0; --i)
            _data[i] = _data[2 * i] + _data[2 * i + 1];
    }

    /// Gets the size of the segment tree.
    [[nodiscard]] constexpr size_t size() const { return _data.size() / 2; }

    /// Gets the data of the segment tree.
    [[nodiscard]] constexpr const std::vector<Node> &data() const { return _data; }

    /// Accesses the element at given position.
    [[nodiscard]] constexpr const Node &operator[](size_t i) const { return _data[i + size()]; }

    /// Updates the element at given position.
    constexpr void update(size_t i, Node node)
    {
        i += size();
        _data[i] = std::move(node);
        for (i /= 2; i; i /= 2)
            _data[i] = _data[2 * i] + _data[2 * i + 1];
    }

    /// Queries the range `[i, j)` of the segment tree.
    [[nodiscard]] constexpr Node query(size_t i, size_t j) const
    {
        i += size(), j += size();
        Node l, r;
        for (; i < j; i /= 2, j /= 2)
        {
            if (i & 1)
                l = l + _data[i++];
            if (j & 1)
                r = _data[--j] + r;
        }
        return std::move(l) + r;
    }

    /// Queries the whole segment tree.
    [[nodiscard]] constexpr Node all() const { return query(0, size()); }

    /// Prints the segment tree.
    friend std::ostream &operator<<(std::ostream &o, const SegmentTree &s)
    {
        o << "Segment tree of size " << s.size() << ":\n  Internal nodes:\n";
        for (size_t i = 1; i < s.size(); ++i)
            o << "    " << i << ": " << s._data[i] << "\n";
        o << "  Leaves:\n";
        for (size_t i = s.size(); i < s._data.size(); ++i)
            o << "    " << i - s.size() << ": " << s._data[i] << "\n";
        return o;
    }
};
} // namespace euler
