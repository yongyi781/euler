#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Unfolds an infinite recursive sequence starting from `a₀` using the function `fn`.
///
/// Usage: `unfold(a₀, fn)` where `fn: T → T` or `fn: T × ℕ → T`.
///
/// For example, `a₁` will be either `fn(a₀)` or `fn(a₀, 1)`.
template <typename T, typename Fn>
    requires(std::invocable<Fn, T> || std::invocable<Fn, T, size_t>)
class unfold : public it_base
{
    T _a0;
    Fn _fn;

  public:
    using value_type = T;

    unfold() = default;
    constexpr unfold(T a0, Fn fn) : _a0(std::move(a0)), _fn(std::move(fn)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        size_t i = 0;
        T state = _a0;
        while (true)
        {
            if (!callbackResult(f, state))
                return result_break;
            if constexpr (std::invocable<Fn, T, size_t>)
            {
                if constexpr (std::same_as<std::invoke_result_t<Fn, T, std::size_t>, std::optional<T>>)
                {
                    auto result = _fn(state, ++i);
                    if (!result)
                        return result_break;
                    state = *result;
                }
                else
                    state = _fn(state, ++i);
            }
            else if constexpr (std::same_as<std::invoke_result_t<Fn, T>, std::optional<T>>)
            {
                auto result = _fn(state);
                if (!result)
                    return result_break;
                state = *result;
            }
            else
            {
                state = _fn(state);
            }
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
};
} // namespace it
} // namespace euler
