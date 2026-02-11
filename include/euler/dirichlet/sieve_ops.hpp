#pragma once

#include "../decls.hpp"

namespace euler::dirichlet
{
namespace detail
{
/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n), using parallelism.
template <typename F, typename G, std::ranges::random_access_range O> void convolveTo_par(F &&f, G &&g, O &&out)
{
    using T = std::ranges::range_value_t<O>;
    if (out.size() <= 1)
        return;
    auto const get_f = [&](size_t i) {
        if constexpr (std::invocable<F, size_t>)
            return f(i);
        else
            return f[i];
    };
    auto const get_g = [&](size_t i) {
        if constexpr (std::invocable<F, size_t>)
            return g(i);
        else
            return g[i];
    };
    size_t const n = out.size() - 1;
    size_t const s = isqrt(n);
    size_t const B = std::max(1UZ, (1UZ << 16) / sizeof(T));
    for (size_t i = 1; i <= s; ++i)
        out[i * i] += T(get_f(i)) * T(get_g(i));
    tbb::parallel_for(
        tbb::blocked_range(1UZ, n + 1, B),
        [&](tbb::blocked_range<size_t> r) {
            size_t const max_i = isqrt(r.end() - 1);
            for (size_t i = 1; i <= max_i; ++i)
            {
                T const a = T(get_f(i)), b = T(get_g(i));
                for (size_t j = std::max(i + 1, (r.begin() + i - 1) / i), k = i * j; k < r.end(); ++j, k += i)
                    out[k] += a * get_g(j) + b * get_f(j);
            }
        },
        tbb::simple_partitioner{});
}

/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n), using parallelism.
template <typename F, typename G> auto convolve_par(F &&f, G &&g, size_t n)
{
    using T = std::remove_cvref_t<std::common_type_t<std::invoke_result_t<F, size_t>, std::invoke_result_t<G, size_t>>>;
    std::vector<T> res(n + 1);
    convolveTo_par(std::forward<F>(f), std::forward<G>(g), res);
    return res;
}

/// Applies `f *= g` in the Dirichlet convolution sense.
template <std::ranges::random_access_range F, typename G> void convolveInPlace_par(F &&f, G &&g)
{
    using T = std::ranges::range_value_t<F>;
    if (f.size() <= 1)
        return;
    auto const get_g = [&](size_t i) {
        if constexpr (std::invocable<G, size_t>)
            return g(i);
        else
            return g[i];
    };
    size_t const n = f.size() - 1;
    size_t const B = std::max(1UZ, (1UZ << 16) / sizeof(T));
    auto const g1 = get_g(1);
    for (size_t end = n; end >= 2; end >>= 1)
    {
        size_t const start = (end >> 1) + 1;
        tbb::parallel_for(
            tbb::blocked_range(start, end + 1, B),
            [&](auto r) {
                size_t const s = isqrt(r.end() - 1);
                if (g1 != 1)
                    for (size_t i = r.begin(); i < r.end(); ++i)
                        f[i] *= g1;
                for (size_t i = r.begin(); i < r.end(); ++i)
                    f[i] += f[1] * get_g(i);
                for (size_t i = 2; i <= s; ++i)
                {
                    T const a = T(f[i]), b = T(get_g(i));
                    for (size_t j = std::max(i + 1, (r.begin() + i - 1) / i), k = i * j; k < r.end(); ++j, k += i)
                        f[k] += a * get_g(j) + b * f[j];
                }
            },
            tbb::simple_partitioner{});
        for (size_t i = isqrt(start - 1) + 1; i * i <= end; ++i)
            f[i * i] += f[i] * get_g(i);
    }
    if (g1 != 1)
        f[1] *= g1;
}

/// Sets `f /= g` in the Dirichlet convolution sense.
template <std::ranges::random_access_range F, typename G> void invConvolve_par(F &&f, G &&g)
{
    using T = std::ranges::range_value_t<F>;
    if (f.size() <= 1)
        return;
    auto const get_g = [&](size_t i) {
        if constexpr (std::invocable<G, size_t>)
            return g(i);
        else
            return g[i];
    };
    size_t const n = f.size() - 1;
    size_t const B = std::max(1UZ, (1UZ << 16) / sizeof(T));
    auto const g1 = get_g(1);
    if (g1 != 1)
        f[1] /= g1;
    for (size_t start = 2; start <= n; start <<= 1)
    {
        size_t const end = std::min(n, (start << 1) - 1);
        size_t const s = isqrt(end);
        for (size_t i = isqrt(start - 1) + 1; i <= s; ++i)
            f[i * i] -= f[i] * get_g(i);
        tbb::parallel_for(
            tbb::blocked_range(start, end + 1, B),
            [&](auto r) {
                for (size_t i = r.begin(); i < r.end(); ++i)
                    f[i] -= f[1] * get_g(i);
                for (size_t i = 2; i <= s; ++i)
                {
                    T const a = T(f[i]), b = T(get_g(i));
                    for (size_t j = std::max(i + 1, (r.begin() + i - 1) / i), k = i * j; k < r.end(); ++j, k += i)
                        f[k] -= a * get_g(j) + b * f[j];
                }
                if (g1 != 1)
                    for (size_t i = r.begin(); i < r.end(); ++i)
                        f[i] /= g1;
            },
            tbb::simple_partitioner{});
    }
}
} // namespace detail

/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n).
template <typename F, typename G, std::ranges::random_access_range O> void convolveInto(F &&f, G &&g, O &&out)
{
    using T = std::ranges::range_value_t<O>;
    if (out.size() <= 1)
        return;
    auto const get_f = [&](size_t i) {
        if constexpr (std::invocable<F, size_t>)
            return f(i);
        else
            return f[i];
    };
    auto const get_g = [&](size_t i) {
        if constexpr (std::invocable<F, size_t>)
            return g(i);
        else
            return g[i];
    };
    size_t const n = out.size() - 1;
    size_t const s = isqrt(n);
    size_t const B = std::max(1UZ, (1UZ << 17) / sizeof(T));
    for (size_t i = 1; i <= s; ++i)
        out[i * i] += T(get_f(i)) * T(get_g(i));
    for (size_t start = 1; start <= n; start += B)
    {
        size_t const end = std::min(n, start + B - 1);
        size_t const max_i = isqrt(end);
        for (size_t i = 1; i <= max_i; ++i)
        {
            T const a = T(get_f(i)), b = T(get_g(i));
            for (size_t j = std::max(i + 1, (start + i - 1) / i), k = i * j; k <= end; ++j, k += i)
                out[k] += a * get_g(j) + b * get_f(j);
        }
    }
}

/// Sets `f *= g` in the Dirichlet convolution sense.
template <std::ranges::random_access_range F, typename G> void convolveInPlace(F &&f, G &&g)
{
    using T = std::ranges::range_value_t<F>;
    if (f.size() <= 1)
        return;
    auto const get_g = [&](size_t i) {
        if constexpr (std::invocable<G, size_t>)
            return g(i);
        else
            return g[i];
    };
    size_t const n = f.size() - 1;
    size_t const B = std::max(1UZ, (1UZ << 17) / sizeof(T));
    auto const g1 = get_g(1);
    for (size_t end = n; end >= 2; end >>= 1)
    {
        size_t const start = (end >> 1) + 1;
        size_t const s = isqrt(end);
        for (size_t min_i = start; min_i <= end; min_i += B)
        {
            size_t const max_i = std::min(end, min_i + B - 1);
            if (g1 != 1)
                for (size_t i = min_i; i <= max_i; ++i)
                    f[i] *= g1;
            for (size_t i = min_i; i <= max_i; ++i)
                f[i] += f[1] * get_g(i);
            for (size_t i = 2; i <= s; ++i)
            {
                T const a = T(f[i]), b = T(get_g(i));
                for (size_t j = std::max(i + 1, (min_i + i - 1) / i), k = i * j; k <= max_i; ++j, k += i)
                    f[k] += a * get_g(j) + b * f[j];
            }
        }
        for (size_t i = isqrt(start - 1) + 1; i <= s; ++i)
            f[i * i] += f[i] * get_g(i);
    }
    if (g1 != 1)
        f[1] *= g1;
}

/// Sets `f /= g` in the Dirichlet convolution sense.
template <std::ranges::random_access_range F, typename G> void invConvolve(F &&f, G &&g)
{
    using T = std::ranges::range_value_t<F>;
    if (f.size() <= 1)
        return;
    auto const get_g = [&](size_t i) {
        if constexpr (std::invocable<G, size_t>)
            return g(i);
        else
            return g[i];
    };
    size_t const n = f.size() - 1;
    size_t const B = std::max(1UZ, (1UZ << 17) / sizeof(T));
    auto const g1 = get_g(1);
    if (g1 != 1)
        f[1] /= g1;
    for (size_t start = 2; start <= n; start <<= 1)
    {
        size_t const end = std::min(n, (start << 1) - 1);
        size_t const s = isqrt(end);
        for (size_t i = isqrt(start - 1) + 1; i <= s; ++i)
            f[i * i] -= f[i] * get_g(i);
        for (size_t min_i = start; min_i <= end; min_i += B)
        {
            size_t const max_i = std::min(end, min_i + B - 1);
            for (size_t i = min_i; i <= max_i; ++i)
                f[i] -= f[1] * get_g(i);
            for (size_t i = 2; i <= s; ++i)
            {
                T const a = T(f[i]), b = T(get_g(i));
                for (size_t j = std::max(i + 1, (min_i + i - 1) / i), k = i * j; k <= max_i; ++j, k += i)
                    f[k] -= a * get_g(j) + b * f[j];
            }
            if (g1 != 1)
                for (size_t i = min_i; i <= max_i; ++i)
                    f[i] /= g1;
        }
    }
}

/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n).
template <execution_policy Exec, typename F, typename G, std::ranges::random_access_range O>
void convolveInto(Exec && /*exec*/, F &&f, G &&g, O &&out)
{
    if constexpr (parallel_policy<Exec>)
        detail::convolveTo_par(std::forward<F>(f), std::forward<G>(g), out);
    else
        convolveInto(std::forward<F>(f), std::forward<G>(g), out);
}

/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n).
template <typename F, typename G> auto convolve(F &&f, G &&g, size_t n)
{
    using T = std::remove_cvref_t<std::common_type_t<std::invoke_result_t<F, size_t>, std::invoke_result_t<G, size_t>>>;
    std::vector<T> res(n + 1);
    convolveInto(std::forward<F>(f), std::forward<G>(g), res);
    return res;
}

/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n), with the specified execution policy.
template <execution_policy Exec, typename F, typename G> auto convolve(Exec && /*exec*/, F &&f, G &&g, size_t n)
{
    if constexpr (parallel_policy<Exec>)
        return detail::convolve_par(std::forward<F>(f), std::forward<G>(g), n);
    else
        return convolve(std::forward<F>(f), std::forward<G>(g), n);
}

/// Computes a sieve for the function `k ↦ ∑ (a*b = k) f(a) * g(b)` in O(n log n), with the specified execution policy.
template <execution_policy Exec, std::ranges::random_access_range F, typename G>
void convolveInPlace(Exec && /*exec*/, F &&f, G &&g)
{
    if constexpr (parallel_policy<Exec>)
        return detail::convolveInPlace_par(f, std::forward<G>(g));
    else
        return convolveInPlace(f, std::forward<G>(g));
}

/// Sets f to f / g in the Dirichlet convolution sense, with the specified execution policy.
template <execution_policy Exec, std::ranges::random_access_range F, typename Fun>
void invConvolve(Exec && /*exec*/, F &&f, Fun g)
{
    if constexpr (parallel_policy<Exec>)
        detail::invConvolve_par(f, std::move(g));
    else
        invConvolve(f, std::move(g));
}
} // namespace euler::dirichlet
