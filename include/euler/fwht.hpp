#pragma once

#include "euler/ZMod.hpp"
#include <ranges>

inline namespace euler
{
/// Perform the Walsh-Hadamard transform in-place.
template <std::ranges::random_access_range Range> void fwht(Range &v)
{
    using T = Range::value_type;
    assert(std::has_single_bit(v.size()));
    for (size_t h = 1; h < v.size(); h *= 2)
        for (size_t i = 0; i < v.size(); i += 2 * h)
            for (size_t j = i; j < i + h; ++j)
                std::tie(v[j], v[j + h]) = std::pair{T(v[j] + v[j + h]), T(v[j] - v[j + h])};
}

/// Perform the inverse Fast Walsh-Hadamard transform in-place.
template <std::ranges::random_access_range Range> void ifwht(Range &v)
{
    assert(std::has_single_bit(v.size()));
    for (size_t h = 1; h < v.size(); h *= 2)
        for (size_t i = 0; i < v.size(); i += 2 * h)
            for (size_t j = i; j < i + h; ++j)
                std::tie(v[j], v[j + h]) = std::pair{(v[j] + v[j + h]) / 2, (v[j] - v[j + h]) / 2};
}

/// Perform the inverse Fast Walsh-Hadamard transform in-place, with an odd modulus.
template <integral2 auto M>
    requires(M % 2 == 1)
void ifwht(std::vector<ZMod<M>> &v)
{
    assert(std::has_single_bit(v.size()));
    auto const inv2 = (M + 1) / 2;
    for (size_t h = 1; h < v.size(); h *= 2)
        for (size_t i = 0; i < v.size(); i += 2 * h)
            for (size_t j = i; j < i + h; ++j)
                std::tie(v[j], v[j + h]) = std::pair{(v[j] + v[j + h]) * inv2, (v[j] - v[j + h]) * inv2};
}

/// Perform the Walsh-Hadamard OR transform in-place.
template <std::ranges::random_access_range Range> void fwht_or(Range &v)
{
    assert(std::has_single_bit(v.size()));
    for (size_t h = 1; h < v.size(); h *= 2)
        for (size_t i = 0; i < v.size(); i += 2 * h)
            for (size_t j = i; j < i + h; ++j)
                v[j + h] += v[j];
}

/// Perform the inverse Fast Walsh-Hadamard OR transform in-place.
template <std::ranges::random_access_range Range> void ifwht_or(Range &v)
{
    assert(std::has_single_bit(v.size()));
    for (size_t h = 1; h < v.size(); h *= 2)
        for (size_t i = 0; i < v.size(); i += 2 * h)
            for (size_t j = i; j < i + h; ++j)
                v[j + h] -= v[j];
}
} // namespace euler
