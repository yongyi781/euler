#pragma once

#include "math.hpp"

namespace euler
{
/// Generate the B + 1 numbers (B(1) = 1/2) up to `k`.
template <typename T> std::vector<T> bernoulliPlus(size_t k)
{
    std::vector<T> A(k + 1), B(k + 1);
    for (size_t m = 0; m <= k; m++)
    {
        A[m] = T(1) / T(m + 1);
        for (int j = m; j >= 1; j--)
            A[j - 1] = j * (A[j - 1] - A[j]);
        if (m == 1 || m % 2 == 0)
            B[m] = A[0];
    }
    return B;
}

/// Stores Bernoulli+ numbers up to k, with the ability to compute power sums up to k.
template <typename T> class Bernoulli
{
    std::vector<T> B;
    std::vector<T> invs;
    std::vector<std::vector<T>> binom_table;

  public:
    using value_type = T;

    Bernoulli(int k)
        : B(bernoulliPlus<T>(k)), invs(range(0, k + 1, [](int i) { return i == 0 ? T{} : T(1) / T(i); })),
          binom_table(euler::binomTable<T>(k + 1))
    {
    }

    /// Accessor for Bernoulli numbers.
    [[nodiscard]] constexpr const value_type &operator[](size_t i) const { return B[i]; }

    /// The list of Bernoulli numbers.
    [[nodiscard]] constexpr const std::vector<value_type> &list() const { return B; }
    /// The list whos ith element is 1/i.
    [[nodiscard]] constexpr const std::vector<value_type> &inverses() const { return invs; }
    /// The precomputed binomial table.
    [[nodiscard]] constexpr const std::vector<std::vector<value_type>> &binomTable() const { return binom_table; }

    /// Returns the size of the Bernoulli list.
    [[nodiscard]] constexpr size_t size() const { return B.size(); }

    /// Sums the kth powers 1^k + ... + n^k using Faulhaber's formula.
    template <typename U> [[nodiscard]] constexpr std::common_type_t<value_type, U> powerSum(const U &n, int k) const
    {
        using Tp = std::common_type_t<value_type, U>;
        assert(std::cmp_less(k, B.size()));
        Tp res = 0;
        Tp x = n;
        for (int j = k; j >= 0; --j, x *= n)
            res += Tp(binom_table[k + 1][j] * B[j]) * x;
        res *= invs[k + 1];
        return res;
    }

    /// Returns the polynomial 1^k + ... + x^k using Faulhaber's formula.
    std::vector<value_type> powerSumPoly(int k) const
    {
        assert(std::cmp_less(k, B.size()));
        std::vector<value_type> res(k + 2);
        for (int j = k; j >= 0; --j)
            res[k + 1 - j] = binom_table[k + 1][j] * B[j] * invs[k + 1];
        return res;
    }

    /// Computes p(1) + p(2) + ... + p(x).
    template <typename Poly> std::vector<value_type> faulhaber(const Poly &p) const
    {
        std::vector<value_type> res(p.size() + 1);
        value_type a;
        for (size_t k = 0; k < p.size(); ++k)
        {
            a = p[k] * invs[k + 1];
            for (size_t j = 0; j <= k; ++j)
                res[k - j + 1] += a * binom_table[k + 1][j] * B[j];
        }
        return res;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                         const Bernoulli<value_type> &b)
    {
        for (size_t i = 0; i < b.B.size(); ++i)
            o << "B[" << i << "] = " << b[i] << "\n";
        return o;
    }
};
} // namespace euler
