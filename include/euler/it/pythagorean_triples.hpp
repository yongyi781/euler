#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates primitive Pythaogrean triples. Takes about 7 ns per triple. Fast fast!
template <integral2 T = int64_t> class primitive_pythagorean_triples : public it_base
{
  public:
    using value_type = std::array<T, 3>;

    primitive_pythagorean_triples() = default;
    constexpr primitive_pythagorean_triples(T limit) : _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return _enumerate(T(3), T(4), T(5), f);
    }

    [[nodiscard]] constexpr T limit() const { return _limit; }

  private:
    T _limit;

    template <std::invocable<value_type> Fun> constexpr result_t _enumerate(T a, T b, T c, Fun f) const
    {
        if (c > _limit)
            return result_continue;
        result_t result = callbackResult(f, value_type{a, b, c});
        if (result == result_break)
            return result_break;
        if (result == result_continue_no_recurse)
            return result_continue;
        std::array arr{std::array{2 * a + b - c, -2 * a + 2 * b + 2 * c, -2 * a + b + 3 * c},
                       std::array{2 * a + b + c, 2 * a - 2 * b + 2 * c, 2 * a - b + 3 * c},
                       std::array{2 * a - b + c, 2 * a + 2 * b + 2 * c, 2 * a + b + 3 * c}};
        for (auto &&[x, y, z] : arr)
        {
            result = _enumerate(x, y, z, f);
            if (result == result_break)
                return result_break;
            if (result == result_continue_no_recurse)
                return result_continue;
        }
        return result_continue;
    }
};

/// Enumerates primitive Pythaogrean triples. Takes about 27 ns per triple.
template <integral2 T = int64_t> class primitive_pythagorean_triples2 : public it_base
{
  public:
    using value_type = std::array<T, 3>;

    primitive_pythagorean_triples2() = default;
    constexpr primitive_pythagorean_triples2(T limit) : _limit(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        using euler::gcd;
        using std::gcd;

        T hv = isqrt(_limit / 2);
        for (T v = 1; v <= hv; ++v)
        {
            for (T u = v + 1; u * u + v * v <= _limit; u += 2)
            {
                if (gcd(u, v) != 1)
                    continue;
                if (!callbackResult(f, value_type{u * u - v * v, 2 * u * v, u * u + v * v}))
                    return result_break;
            }
        }
        return result_continue;
    }

    constexpr T limit() const { return _limit; }

  private:
    T _limit;
};

/// Enumerates Pythaogrean triples. Takes about 0.5 ns per triple if using the faster version
/// of the primitive triple generator. Fast fast!
template <typename T = int64_t, typename P = primitive_pythagorean_triples<T>>
class pythagorean_triples : public it_base
{
  public:
    using value_type = P::value_type;

    pythagorean_triples() = default;
    constexpr pythagorean_triples(T limit) : base(limit) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        return base([&](auto &&x) {
            auto &&[a, b, c] = std::forward<decltype(x)>(x);
            T hk = base.limit() / c;
            for (T k = 1; k <= hk; ++k)
                if (!callbackResult(f, value_type{k * a, k * b, k * c}))
                    return result_break;
            return result_continue;
        });
    }

  private:
    P base;
};
} // namespace it
} // namespace euler
