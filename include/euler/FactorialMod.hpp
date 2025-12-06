#pragma once

#include "ZMod.hpp"

namespace euler
{
/// Class to compute factorials and binomial coefficients modulo a compile-time constant modulus.
/// Modulus is a template parameter for compile time optimization benefits.
template <integral2 auto M, bool ComputeInverse = true> class FactorialMod
{
    using T = decltype(M);

    std::vector<ZMod<M>> fact;
    std::vector<ZMod<M>> invFact;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        fact[0] = 1;
        std::transform_inclusive_scan(std::forward<Exec>(exec), counting_iterator(T(1)),
                                      counting_iterator(T(fact.size())), std::next(fact.begin()), std::multiplies{},
                                      [](T x) { return ZMod<M>(x); });
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(std::forward<Exec>(exec), counting_iterator(T(0)),
                                      counting_iterator(T(invFact.size())), invFact.rbegin(), ~fact.back(),
                                      std::multiplies{}, [&](T x) { return ZMod<M>(invFact.size() - 1 - x); });
    }

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
};

/// Class to compute factorials and binomial coefficients modulo a dynamic modulus.
template <integral2 T, bool ComputeInverse = true> class FactorialModX
{
    T p_{};
    std::vector<T> fact;
    std::vector<T> invFact;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        fact[0] = 1;
        std::inclusive_scan(std::forward<Exec>(exec), counting_iterator(T(1)), counting_iterator(T(fact.size())),
                            std::next(fact.begin()), mod_multiplies{p_});
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(
            std::forward<Exec>(exec), counting_iterator(T(0)), counting_iterator(T(invFact.size())), invFact.rbegin(),
            modInverse(fact.back(), p_), mod_multiplies{p_}, [&](T x) -> T { return invFact.size() - 1 - x; });
    }

  public:
    using value_type = T;

    template <execution_policy Exec>
        requires(!std::is_integral_v<Exec>)
    FactorialModX(Exec &&exec, T p, size_t limit) : p_(p), fact(limit + 1, 0), invFact(limit + 1, 0)
    {
        if (limit > 0)
        {
            computeFactorials(std::forward<Exec>(exec));
            if constexpr (ComputeInverse)
                computeInverseFactorials(std::forward<Exec>(exec));
        }
    }

    explicit FactorialModX(T p, size_t limit) : FactorialModX(std::execution::seq, p, limit) {}

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
            return modInverse(fact[n], p_);
    }

    [[nodiscard]] value_type binomial(int64_t n, int64_t r) const
    {
        assert(n >= 0 && n <= (int64_t)fact.size());
        if (r < 0 || r > n)
            return (int64_t)0;
        return factorial(n) * inverseFactorial(n - r) % p_ * inverseFactorial(r) % p_;
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
};
} // namespace euler
