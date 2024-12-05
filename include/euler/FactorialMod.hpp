#pragma once

#include "ZMod.hpp"

inline namespace euler
{
/// Class to compute factorials and binomial coefficients modulo a compile-time constant modulus.
/// Modulus is a template parameter for compile time optimization benefits.
template <integral2 auto M, bool ComputeInverse = true> class FactorialMod
{
  public:
    using value_type = ZMod<M>;

    template <execution_policy Exec>
        requires(!std::is_integral_v<Exec>)
    FactorialMod(Exec &&exec, size_t limit) : fact(limit + 1, 0), invFact(limit + 1, 0)
    {
        if (limit > 0)
        {
            computeFactorials(std::forward<Exec>(exec));
            if constexpr (ComputeInverse)
                computeInverseFactorials(std::forward<Exec>(exec));
        }
    }

    explicit FactorialMod(size_t limit = 0) : FactorialMod(std::execution::seq, limit) {}

    [[nodiscard]] value_type operator()(int64_t n) const { return factorial(n); }

    [[nodiscard]] value_type factorial(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)fact.size());
        return fact[n];
    }

    [[nodiscard]] value_type inverseFactorial(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)fact.size());
        if constexpr (ComputeInverse)
            return invFact[n];
        else
            return ~fact[n];
    }

    [[nodiscard]] value_type binomial(int64_t n, int64_t r) const
    {
        assert(n >= 0 && n <= (int64_t)fact.size());
        if (r < 0 || r > n)
            return (int64_t)0;
        return factorial(n) * inverseFactorial(n - r) * inverseFactorial(r);
    }

    template <execution_policy Exec> void resize(Exec &&exec, size_t limit)
    {
        if (limit > fact.size())
        {
            fact.resize(limit + 1);
            invFact.resize(limit + 1);
            computeFactorials(std::forward<Exec>(exec));
            if constexpr (ComputeInverse)
                computeInverseFactorials(std::forward<Exec>(exec));
        }
    }

    void resize(size_t newSize) { resize(std::execution::par, newSize); }

  private:
    using T = decltype(M);

    std::vector<value_type> fact;
    std::vector<value_type> invFact;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        fact[0] = 1;
        std::transform_inclusive_scan(std::forward<Exec>(exec), counting_iterator(T(1)),
                                      counting_iterator(T(fact.size())), std::next(fact.begin()), std::multiplies{},
                                      [](T x) { return value_type(x); });
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(std::forward<Exec>(exec), counting_iterator(T(0)),
                                      counting_iterator(T(invFact.size())), invFact.rbegin(), ~fact.back(),
                                      std::multiplies{}, [&](T x) { return value_type(invFact.size() - 1 - x); });
    }
};
} // namespace euler
