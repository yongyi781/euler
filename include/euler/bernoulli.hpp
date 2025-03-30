#pragma once

#include "math.hpp"

inline namespace euler
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
  public:
    Bernoulli(int k)
        : B(bernoulliPlus<T>(k)), invs(range(0, k + 1, [](int i) { return i == 0 ? T{} : T(1) / T(i); })),
          btable(binomTable<T>(k + 1))
    {
    }

    /// Accessor for Bernoulli numbers.
    [[nodiscard]] constexpr const T &operator[](size_t i) const { return B[i]; }

    /// Sums the kth powers 1^k + ... + n^k using Faulhaber's formula.
    template <typename U> [[nodiscard]] constexpr std::common_type_t<T, U> powerSum(const U &n, int k) const
    {
        using Tp = std::common_type_t<T, U>;
        assert(std::cmp_less(k, B.size()));
        Tp res = 0;
        Tp x = n;
        for (int j = k; j >= 0; --j, x *= n)
            res += Tp(btable[k + 1][j] * B[j]) * x;
        res *= invs[k + 1];
        return res;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const Bernoulli<T> &b)
    {
        for (size_t i = 0; i < b.B.size(); ++i)
            o << "B[" << i << "] = " << b[i] << "\n";
        return o;
    }

  private:
    std::vector<T> B;
    std::vector<T> invs;
    std::vector<std::vector<T>> btable;
};
} // namespace euler
