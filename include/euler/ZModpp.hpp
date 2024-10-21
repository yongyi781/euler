#pragma once

#include "ZMod.hpp"
#include "prime.hpp"

inline namespace euler
{
/// Specialized class for integers modulo prime powers.
template <std::integral auto P, std::integral auto E>
    requires(isPrime(P))
class ZModpp : public ZMod<pow(P, E)>
{
  public:
    using value_type = decltype(P);

    ZModpp() = default;
    constexpr ZModpp(const value_type &value) : ZMod<pow(P, E)>(value) {}
    constexpr ZModpp(const ZMod<pow(P, E)> &value) : ZMod<pow(P, E)>(value) {}

    constexpr std::optional<ZModpp> sqrt() const
    {
        auto a = this->value();
        if (a == 0)
            return 0;
        auto v = valuationDivide(a, P);
        if (v % 2 != 0)
            return std::nullopt;
        auto x = ZMod<P>(a).sqrt();
        if (!x)
            return std::nullopt;
        auto y = ZModpp(x->value());
        constexpr value_type r = ZModpp::modulus / P;
        auto s = y.pow(r) * ZModpp(a).pow((ZModpp::modulus - 2 * r + 1) / 2) * pow(P, v / 2);
        return s.value() <= ZModpp::modulus / 2 ? s : -s;
    }
};
} // namespace euler
