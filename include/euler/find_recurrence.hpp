#pragma once

#include "algorithm.hpp"
#include "modular_arithmetic.hpp"


/// Finds a linear recurrence for `v` using the Berlekamp-Massey algorithm.
template <integral2 T> std::vector<T> findRecurrence(const std::vector<T> &v, T modulus)
{
    size_t n = v.size();
    size_t L = 0;
    size_t m = 0;
    std::vector<T> C(n);
    std::vector<T> B(n);
    C[0] = B[0] = 1;
    T b = 1;
    for (size_t i = 0; i < n; ++i)
    {
        ++m;
        T d = v[i] % modulus;
        for (size_t j = 1; j <= L; ++j)
            d = (d + C[j] * v[i - j]) % modulus;
        if (d == 0)
            continue;
        auto C2 = C;
        T const coeff = d * modInverse(b, modulus) % modulus;
        for (size_t j = m; j < n; ++j)
            C[j] = (C[j] - coeff * B[j - m]) % modulus;
        if (i < 2 * L)
            continue;
        L = i + 1 - L;
        B = C2;
        b = d;
        m = 0;
    }
    C.resize(L + 1);
    C.erase(C.begin());
    for (auto &x : C)
    {
        x = (modulus - x) % modulus;
        // For the user's benefit
        if (x > modulus / 2 + 1)
            x -= modulus;
    }
    return C;
}

/// Finds a linear recurrence for `v` using the Berlekamp-Massey algorithm. Returns a list `[c1, ..., ck]` such that
/// the charpoly of the recurrence is `x^k - c1 * x^(k-1) + ... + ck`.
/// Note: It is not recommended to call this function with vectors of integral types.
template <std::ranges::range Range> std::vector<std::ranges::range_value_t<Range>> findRecurrence(Range &&v)
{
    using T = std::ranges::range_value_t<Range>;

    size_t const n = v.size();
    size_t L = 0;
    size_t m = 0;
    std::vector<T> C(n);
    std::vector<T> B(n);
    C[0] = B[0] = 1;
    T b = 1;
    for (size_t i = 0; i < n; ++i)
    {
        ++m;
        T const d = sum(0, L, [&](size_t j) { return C[j] * v[i - j]; });
        if (d == 0)
            continue;
        auto C2 = C;
        T const coeff = d / b;
        for (size_t j = m; j < n; ++j)
            C[j] -= coeff * B[j - m];
        if (i < 2 * L)
            continue;
        L = i + 1 - L;
        B = C2;
        b = d;
        m = 0;
    }
    C.resize(L + 1);
    C.erase(C.begin());
    for (T &x : C)
        x = -x;
    return C;
}
