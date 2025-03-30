#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// D-finite recurrence.
template <typename T = int64_t> class dfinite : public it_base
{
  public:
    using value_type = T;

    constexpr dfinite(std::vector<std::vector<int64_t>> coeffPolys, std::vector<T> init, int64_t firstIndex = 0)
        : _coeffPolys(std::move(coeffPolys)), _init(std::move(init)), _firstIndex(firstIndex)
    {
        assert(_coeffPolys.size() - 1 == _init.size() &&
               "There must be the same number of coefficients and initial values.");
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        for (auto &&c : _init)
            if (!callbackResult(f, c))
                return result_break;
        int64_t const order = _coeffPolys.size() - 1;
        std::vector<T> y = _init;
        std::vector<T> z(order);
        for (int64_t n = _firstIndex + order;; ++n)
        {
            std::copy(y.begin() + 1, y.end(), z.begin());
            z.back() = 0;
            for (size_t k = 0; k < y.size(); ++k)
                z.back() -= evalPoly(_coeffPolys[k], n - order) * y[k];
            z.back() /= evalPoly(_coeffPolys.back(), n - order);
            if (!callbackResult(f, z.back()))
                return result_break;
            swap(y, z);
        }
        return result_continue;
    }

    [[nodiscard]] constexpr size_t order() const { return _coeffPolys.size() - 1; }
    [[nodiscard]] constexpr size_t degree() const
    {
        size_t res = 0;
        for (auto &&poly : _coeffPolys)
            res = std::max(res, poly.size());
        return res;
    }

    /// It makes more sense to index it this way.
    [[nodiscard]] constexpr value_type operator[](int64_t i) const
    {
        int64_t n = i - _firstIndex;
        assert(n >= 0 && "Index out of range.");
        if (n < (int64_t)_init.size())
            return _init[n];
        value_type result{};
        (*this)([&](auto &&x) {
            if (n == 0)
            {
                result = std::forward<decltype(x)>(x);
                return result_break;
            }
            --n;
            return result_continue;
        });
        return result;
    }

    /// Gets the ith term of this enumerable, wrapped in a `std::optional`.
    [[nodiscard]] constexpr std::optional<value_type> at(int64_t i) const { return (*this)[i]; }

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr size_t size() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr auto cycle() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr auto last() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr auto reduce() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr auto sum() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr auto product() const = delete;

  private:
    std::vector<std::vector<int64_t>> _coeffPolys;
    std::vector<T> _init;
    int64_t _firstIndex;

    [[nodiscard]] constexpr T evalPoly(const std::vector<int64_t> &poly, const T &x) const
    {
        T res = 0;
        for (auto &&c : poly | std::views::reverse)
        {
            res *= x;
            res += c;
        }
        return res;
    }
};
} // namespace it
} // namespace euler
