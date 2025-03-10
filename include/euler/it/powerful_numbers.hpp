#pragma once

#include "euler/prime.hpp"
#include "tree.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates powerful numbers up to a given N, along with their factorizations.
/// Usage: `powerful_numbers_factored(<N>)`.
class powerful_numbers_factored : public it_base
{
  public:
    using value_type = std::pair<int64_t, Factorization<int64_t>>;
    using value_reference_type = std::pair<int64_t, const Factorization<int64_t> &>;

    powerful_numbers_factored() = default;
    constexpr powerful_numbers_factored(int64_t limit) : _limit(limit), _primes(primeRange(isqrt(limit))) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        Factorization<int64_t> fac{};
        fac.reserve(4);
        return it::tree_preorder(
            std::tuple<int64_t, Factorization<int64_t> &, It>{1_i64, fac, _primes.begin()},
            [&](auto &&x, auto f) {
                auto &&[n, fac, i] = x;
                bool currPrimeInFac = !fac.empty();
                for (auto j = i; j != std::ranges::end(_primes); ++j)
                {
                    int64_t const p = currPrimeInFac ? *j : *j * *j;
                    if (!mulLeq(n, p, _limit))
                        return result_continue;
                    if (currPrimeInFac)
                        ++fac.back().second;
                    else
                        fac.emplace_back(*j, 2);
                    if (!f({n * p, fac, j}))
                        return result_break;
                    if (currPrimeInFac)
                    {
                        --fac.back().second;
                        currPrimeInFac = false;
                    }
                    else
                        fac.pop_back();
                }
                return result_continue;
            },
            [&](auto &&x) {
                auto &&[n, fac, i] = x;
                if (fac.empty() || fac.back().first != *i)
                    return mulLeq(n, *i * *i, _limit);
                return mulLeq(n, *i, _limit);
            })([&](auto &&x) { return f(value_reference_type{get<0>(x), get<1>(x)}); });
    }

  private:
    using It = std::vector<int64_t>::const_iterator;

    int64_t _limit = 0;
    std::vector<int64_t> _primes;
};
} // namespace it
} // namespace euler
