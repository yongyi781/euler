#pragma once

#include "../PF.hpp"
#include "../prime.hpp"
#include "tree.hpp"

namespace euler::it
{
/// Enumerates powerful numbers up to a given N, along with their factorizations.
/// Usage: `powerful_numbers_factored(<N>)`.
class powerful_numbers_factored : public it_base
{
    using It = std::vector<u64>::const_iterator;

    u64 limit_ = 0;
    std::vector<u64> primes_;

  public:
    using value_type = std::pair<u64, PF<u64>>;
    using value_reference_type = std::pair<u64, const PF<u64> &>;

    powerful_numbers_factored() = default;
    constexpr powerful_numbers_factored(u64 limit) : limit_(limit), primes_(primeRange(isqrt(limit))) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        PF<u64> fac{};
        fac.data().reserve(4);
        return it::tree_preorder(
            std::tuple<u64, PF<u64> &, It>{1_u64, fac, primes_.begin()},
            [&](auto &&x, auto f) {
                auto &&[n, fac, i] = x;
                bool currPrimeInFac = !fac.empty();
                for (auto j = i; j != std::ranges::end(primes_); ++j)
                {
                    u64 const p = currPrimeInFac ? *j : *j * *j;
                    if (!mulLeq(n, p, limit_))
                        return result_continue;
                    if (currPrimeInFac)
                        ++fac.back().second;
                    else
                        fac.data().emplace_back(*j, 2);
                    if (!f({n * p, fac, j}))
                        return result_break;
                    if (currPrimeInFac)
                    {
                        --fac.back().second;
                        currPrimeInFac = false;
                    }
                    else
                        fac.data().pop_back();
                }
                return result_continue;
            },
            [&](auto &&x) {
                auto &&[n, fac, i] = x;
                if (fac.empty() || fac.back().first != *i)
                    return mulLeq(n, *i * *i, limit_);
                return mulLeq(n, *i, limit_);
            })([&](auto &&x) { return f(value_reference_type{get<0>(x), get<1>(x)}); });
    }
};
} // namespace euler::it
