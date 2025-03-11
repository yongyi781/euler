#pragma once

#include "math.hpp"

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
    [[nodiscard]] constexpr T powerSum(const T &n, int k) const
    {
        assert(std::cmp_less(k, B.size()));
        T res = 0;
        T x = n;
        for (int j = k; j >= 0; --j, x *= n)
            res += btable[k + 1][j] * B[j] * x;
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
