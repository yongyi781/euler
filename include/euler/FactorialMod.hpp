#pragma once

#include "ZMod.hpp"

namespace euler
{
/// Class to compute factorials and binomial coefficients modulo a compile-time constant modulus.
/// Modulus is a template parameter for compile time optimization benefits.
template <integral2 auto M, bool ComputeInverse = true> class FactorialMod
{
    using T = decltype(M);

    std::vector<ZMod<M>> fact_;
    std::vector<ZMod<M>> ifact_;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        fact_[0] = 1;
        std::transform_inclusive_scan(std::forward<Exec>(exec), counting_iterator(1UZ), counting_iterator(fact_.size()),
                                      std::next(fact_.begin()), std::multiplies{}, [](size_t i) { return ZMod<M>(i); });
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(std::forward<Exec>(exec), counting_iterator(0UZ),
                                      counting_iterator(ifact_.size()), ifact_.rbegin(), ~fact_.back(),
                                      std::multiplies{}, [&](size_t i) { return ZMod<M>(ifact_.size() - i - 1); });
    }

  public:
    using value_type = ZMod<M>;

    template <execution_policy Exec>
        requires(!integral2<Exec>)
    FactorialMod(Exec &&exec, size_t limit) : fact_(limit + 1, 0), ifact_(limit + 1, 0)
    {
        if (limit == 0)
            return;
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    explicit FactorialMod(size_t limit = 0) : FactorialMod(std::execution::seq, limit) {}

    /// Returns n!.
    [[nodiscard]] value_type operator()(int64_t n) const { return get(n); }

    /// Returns n!.
    [[nodiscard]] value_type get(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)fact_.size());
        return fact_[n];
    }

    /// Returns 1 / n!.
    [[nodiscard]] value_type inv(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)fact_.size());
        if constexpr (ComputeInverse)
            return ifact_[n];
        else
            return ~fact_[n];
    }

    [[nodiscard]] value_type binomial(int64_t n, int64_t r) const
    {
        assert(n >= 0 && n <= (int64_t)fact_.size());
        if (r < 0 || r > n)
            return (int64_t)0;
        if constexpr (ComputeInverse)
            return fact_[n] * ifact_[n - r] * ifact_[r];
        else
            return fact_[n] / (fact_[n - r] * fact_[r]);
    }

    template <execution_policy Exec> void resize(Exec &&exec, size_t limit)
    {
        if (limit <= fact_.size())
            return;
        fact_.resize(limit + 1);
        ifact_.resize(limit + 1);
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    void resize(size_t newSize) { resize(std::execution::par, newSize); }

    const std::vector<value_type> &fact() const { return fact_; }
    const std::vector<value_type> &ifact() const { return ifact_; }
};

/// Class to compute factorials and binomial coefficients modulo a dynamic modulus.
template <integral2 T, bool ComputeInverse = true> class FactorialModX
{
    T p_{};
    std::vector<T> fact_;
    std::vector<T> ifact_;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        fact_[0] = 1;
        std::inclusive_scan(std::forward<Exec>(exec), counting_iterator(T(1)), counting_iterator(T(fact_.size())),
                            std::next(fact_.begin()), mod_multiplies{p_});
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(
            std::forward<Exec>(exec), counting_iterator(T(0)), counting_iterator(T(ifact_.size())), ifact_.rbegin(),
            modInverse(fact_.back(), p_), mod_multiplies{p_}, [&](T x) -> T { return ifact_.size() - x - 1; });
    }

  public:
    using value_type = T;

    template <execution_policy Exec>
        requires(!integral2<Exec>)
    FactorialModX(Exec &&exec, T p, size_t limit) : p_(p), fact_(limit + 1, 0), ifact_(limit + 1, 0)
    {
        if (limit == 0)
            return;
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    explicit FactorialModX(T p, size_t limit = 0) : FactorialModX(std::execution::seq, p, limit) {}

    /// Returns n!.
    [[nodiscard]] value_type operator()(int64_t n) const { return get(n); }

    /// Returns n!.
    [[nodiscard]] value_type get(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)fact_.size());
        return fact_[n];
    }

    /// Returns 1 / n!.
    [[nodiscard]] value_type inv(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)fact_.size());
        if constexpr (ComputeInverse)
            return ifact_[n];
        else
            return modInverse(fact_[n], p_);
    }

    [[nodiscard]] value_type binomial(int64_t n, int64_t r) const
    {
        assert(n >= 0 && n <= (int64_t)fact_.size());
        if (r < 0 || r > n)
            return 0;
        return fact_[n] * ifact_[n - r] % p_ * ifact_[r] % p_;
    }

    template <execution_policy Exec> void resize(Exec &&exec, size_t limit)
    {
        if (limit <= fact_.size())
            return;
        fact_.resize(limit + 1);
        ifact_.resize(limit + 1);
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    void resize(size_t newSize) { resize(std::execution::par, newSize); }

    const std::vector<value_type> &fact() const { return fact_; }
    const std::vector<value_type> &ifact() const { return ifact_; }
};
} // namespace euler
