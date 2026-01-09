#pragma once

#include "ZMod.hpp"

namespace euler
{
/// Class to compute factorials and binomial coefficients modulo a compile-time constant modulus.
/// Modulus is a template parameter for compile time optimization benefits.
template <integral2 auto M, bool ComputeInverse = true> class FactorialMod
{
    using T = decltype(M);

    std::vector<ZMod<M>> facts_;
    std::vector<ZMod<M>> ifacts_;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        facts_[0] = 1;
        std::transform_inclusive_scan(std::forward<Exec>(exec), counting_iterator(1UZ),
                                      counting_iterator(facts_.size()), std::next(facts_.begin()), std::multiplies{},
                                      [](size_t i) { return ZMod<M>(i); });
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(std::forward<Exec>(exec), counting_iterator(0UZ),
                                      counting_iterator(ifacts_.size()), ifacts_.rbegin(), ~facts_.back(),
                                      std::multiplies{}, [&](size_t i) { return ZMod<M>(ifacts_.size() - i - 1); });
    }

  public:
    using value_type = ZMod<M>;

    template <execution_policy Exec>
        requires(!integral2<Exec>)
    FactorialMod(Exec &&exec, size_t limit) : facts_(limit + 1, 0), ifacts_(limit + 1, 0)
    {
        if (limit == 0)
            return;
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    explicit FactorialMod(size_t limit = 0) : FactorialMod(std::execution::seq, limit) {}

    /// Returns `n!`.
    [[nodiscard]] value_type operator()(int64_t n) const { return fact(n); }

    /// Returns `n!`.
    [[nodiscard]] value_type fact(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)facts_.size());
        return facts_[n];
    }

    /// Returns `1 / n!`.
    [[nodiscard]] value_type ifact(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)facts_.size());
        if constexpr (ComputeInverse)
            return ifacts_[n];
        else
            return ~facts_[n];
    }

    [[nodiscard]] value_type binomial(int64_t n, int64_t r) const
    {
        assert(n >= 0 && n <= (int64_t)facts_.size());
        if (r < 0 || r > n)
            return (int64_t)0;
        if constexpr (ComputeInverse)
            return facts_[n] * ifacts_[n - r] * ifacts_[r];
        else
            return facts_[n] / (facts_[n - r] * facts_[r]);
    }

    /// Returns `1 / n` (mod `p`).
    [[nodiscard]] value_type inv(int64_t n) const { return ifacts_[n] * facts_[n - 1]; }

    template <execution_policy Exec> void resize(Exec &&exec, size_t limit)
    {
        if (limit <= facts_.size())
            return;
        facts_.resize(limit + 1);
        ifacts_.resize(limit + 1);
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    void resize(size_t newSize) { resize(std::execution::seq, newSize); }

    const std::vector<value_type> &facts() const { return facts_; }
    const std::vector<value_type> &ifacts() const { return ifacts_; }
};

/// Class to compute factorials and binomial coefficients modulo a dynamic modulus.
template <integral2 T, bool ComputeInverse = true> class FactorialModX
{
    T p_{};
    std::vector<T> facts_;
    std::vector<T> ifacts_;

    template <execution_policy Exec> void computeFactorials(Exec &&exec)
    {
        facts_[0] = 1;
        std::inclusive_scan(std::forward<Exec>(exec), counting_iterator(T(1)), counting_iterator(T(facts_.size())),
                            std::next(facts_.begin()), mod_multiplies{p_});
    }

    template <execution_policy Exec> void computeInverseFactorials(Exec &&exec)
    {
        std::transform_exclusive_scan(
            std::forward<Exec>(exec), counting_iterator(T(0)), counting_iterator(T(ifacts_.size())), ifacts_.rbegin(),
            modInverse(facts_.back(), p_), mod_multiplies{p_}, [&](T x) -> T { return ifacts_.size() - x - 1; });
    }

  public:
    using value_type = T;

    template <execution_policy Exec>
        requires(!integral2<Exec>)
    FactorialModX(Exec &&exec, T p, size_t limit) : p_(p), facts_(limit + 1, 0), ifacts_(limit + 1, 0)
    {
        if (limit == 0)
            return;
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    explicit FactorialModX(T p, size_t limit = 0) : FactorialModX(std::execution::seq, p, limit) {}

    /// Returns n!.
    [[nodiscard]] value_type operator()(int64_t n) const { return fact(n); }

    /// Returns n!.
    [[nodiscard]] value_type fact(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)facts_.size());
        return facts_[n];
    }

    /// Returns 1 / n!.
    [[nodiscard]] value_type ifact(int64_t n) const
    {
        assert(n >= 0 && n <= (int64_t)facts_.size());
        if constexpr (ComputeInverse)
            return ifacts_[n];
        else
            return modInverse(facts_[n], p_);
    }

    [[nodiscard]] value_type binomial(int64_t n, int64_t r) const
    {
        assert(n >= 0 && n <= (int64_t)facts_.size());
        if (r < 0 || r > n)
            return 0;
        return facts_[n] * ifacts_[n - r] % p_ * ifacts_[r] % p_;
    }

    /// Returns `1 / n` (mod `p`).
    [[nodiscard]] value_type inv(int64_t n) const { return ifacts_[n] * facts_[n - 1] % p_; }

    template <execution_policy Exec> void resize(Exec &&exec, size_t limit)
    {
        if (limit <= facts_.size())
            return;
        facts_.resize(limit + 1);
        ifacts_.resize(limit + 1);
        computeFactorials(std::forward<Exec>(exec));
        if constexpr (ComputeInverse)
            computeInverseFactorials(std::forward<Exec>(exec));
    }

    void resize(size_t newSize) { resize(std::execution::seq, newSize); }

    const std::vector<value_type> &facts() const { return facts_; }
    const std::vector<value_type> &ifacts() const { return ifacts_; }
};
} // namespace euler
