#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// The Fibonacci sequence with initial values `a` (default: 0) and `b` (default: 1).
template <typename T = int64_t> class fibonacci : public it_base
{
  public:
    using value_type = T;

    constexpr fibonacci(T a = T(0), T b = T(1)) : _a(a), _b(b) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        T x = _a;
        T y = _b;

        while (true)
        {
            if (!callbackResult(f, x))
                return result_break;
            std::tie(x, y) = std::pair{y, T(x + y)};
        }
    }

    /// Deleted because this enumerable is guaranteed to be infinite.
    [[nodiscard]] constexpr size_t size() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto cycle() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto last() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto reduce() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto sum() const = delete;

    /// Deleted because this enumerable is guaranteed to be infinite.
    constexpr auto product() const = delete;

  private:
    T _a, _b;
};
} // namespace it
} // namespace euler
