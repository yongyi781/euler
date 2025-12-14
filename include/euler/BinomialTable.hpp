#pragma once

#include "types.hpp"

namespace euler
{
/// Binomial table.
template <typename T = i64> class BinomialTable
{
    size_t max_n_ = 0;
    std::vector<T> data_{1};

    [[nodiscard]] constexpr size_t index(size_t n, size_t k) const { return (n * (n + 1)) / 2 + k; }

  public:
    BinomialTable() = default;
    constexpr explicit BinomialTable(size_t max_n) : max_n_(max_n), data_((max_n + 1) * (max_n + 2) / 2)
    {
        for (size_t n = 0; n <= max_n; ++n)
        {
            data_[index(n, 0)] = 1;
            data_[index(n, n)] = 1;
            for (size_t k = 1; k < n; ++k)
                data_[index(n, k)] = data_[index(n - 1, k - 1)] + data_[index(n - 1, k)];
        }
    }

    [[nodiscard]] constexpr size_t max_n() const { return max_n_; }
    [[nodiscard]] constexpr const std::vector<T> &data() const { return data_; }
    [[nodiscard]] constexpr size_t size() const { return data_.size(); }
    [[nodiscard]] constexpr auto begin() const { return data_.begin(); }
    [[nodiscard]] constexpr auto end() const { return data_.end(); }

    /// Accessor with bounds checking.
    [[nodiscard]] constexpr T operator()(size_t n, size_t k) const
    {
        assert(n <= max_n_);
        if (k > n)
            return 0;
        return data_[index(n, k)];
    }

    /// Returns a span representing the n-th row: [nC0, nC1, ... nCn].
    [[nodiscard]] std::span<const T> operator[](size_t n) const
    {
        assert(n <= max_n_);
        return std::span<const T>(data_.data() + n * (n + 1) / 2, n + 1);
    }

    /// Resizes the table.
    void resize(size_t new_max_n)
    {
        if (new_max_n == max_n_)
            return;
        size_t new_size = (new_max_n + 1) * (new_max_n + 2) / 2;
        if (new_max_n < max_n_)
        {
            data_.resize(new_size);
            max_n_ = new_max_n;
            return;
        }
        data_.resize(new_size);
        for (size_t n = max_n_ + 1; n <= new_max_n; ++n)
        {
            data_[index(n, 0)] = 1;
            data_[index(n, n)] = 1;
            for (size_t k = 1; k < n; ++k)
                data_[index(n, k)] = data_[index(n - 1, k - 1)] + data_[index(n - 1, k)];
        }
        max_n_ = new_max_n;
    }
};
} // namespace euler
