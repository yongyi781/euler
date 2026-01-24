#pragma once

#include <ostream>
#include <ranges>
#include <vector>

namespace euler
{
/// A segment tree for an arbitrary monoid.
template <typename Node, typename BinaryOp = std::plus<>> class SegmentTree
{
    std::vector<Node> data_;
    BinaryOp op_;

  public:
    /// Creates an empty segment tree of given size.
    constexpr explicit SegmentTree(size_t n, BinaryOp op = {}) : data_(2 * n), op_(std::move(op)) {}

    /// Creates an empty segment tree of given size.
    constexpr SegmentTree(size_t n, const Node &value, BinaryOp op = {}) : SegmentTree(n, std::move(op))
    {
        std::ranges::fill(data_.begin() + size(), data_.end(), value);
        for (size_t i = size() - 1; i > 0; --i)
            data_[i] = op_(data_[2 * i], data_[2 * i + 1]);
    }

    /// Creates a segment tree from a range of elements.
    template <std::ranges::sized_range Range>
    constexpr SegmentTree(Range &&init, BinaryOp op = {}) : SegmentTree(init.size(), std::move(op))
    {
        std::ranges::move(std::forward<Range>(init), data_.begin() + size());
        for (size_t i = size() - 1; i > 0; --i)
            data_[i] = op_(data_[2 * i], data_[2 * i + 1]);
    }

    /// Gets the size of the segment tree.
    [[nodiscard]] constexpr size_t size() const { return data_.size() / 2; }

    /// Gets the data of the segment tree.
    [[nodiscard]] constexpr const std::vector<Node> &data() const { return data_; }

    /// Accesses the element at given position.
    [[nodiscard]] constexpr const Node &operator[](size_t i) const { return data_[i + size()]; }

    /// Updates the element at given position.
    constexpr void update(size_t i, Node node)
    {
        i += size();
        data_[i] = std::move(node);
        for (i >>= 1; i; i >>= 1)
            data_[i] = op_(data_[2 * i], data_[2 * i + 1]);
    }

    /// Queries the range `[i, j)` of the segment tree.
    [[nodiscard]] constexpr Node query(size_t i, size_t j) const
    {
        i += size(), j += size();
        Node l, r;
        for (; i < j; i >>= 1, j >>= 1)
        {
            if (i & 1)
                l = op_(l, data_[i++]);
            if (j & 1)
                r = op_(data_[--j], r);
        }
        return op_(std::move(l), std::move(r));
    }

    /// Queries the whole segment tree.
    [[nodiscard]] constexpr Node all() const { return query(0, size()); }

    /// Prints the segment tree.
    friend std::ostream &operator<<(std::ostream &o, const SegmentTree &s)
    {
        o << "Segment tree of size " << s.size() << ":\n  Internal nodes:\n";
        for (size_t i = 1; i < s.size(); ++i)
            o << "    " << i << ": " << s.data_[i] << "\n";
        o << "  Leaves:\n";
        for (size_t i = s.size(); i < s.data_.size(); ++i)
            o << "    " << i - s.size() << ": " << s.data_[i] << "\n";
        return o;
    }
};
} // namespace euler
