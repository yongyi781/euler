#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Enumerates all product partitions of `n`.
template <integral2 T> class product_partitions : public it_base
{
  public:
    using value_type = std::vector<T>;

    product_partitions() = default;
    constexpr product_partitions(T n) : _n(std::move(n)) {}

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        value_type partition;
        return _enumerate(partition, _n, f);
    }

  private:
    T _n;

    template <std::invocable<value_type> Fun> constexpr result_t _enumerate(value_type &partition, T n, Fun f) const
    {
        partition.push_back(n);
        if (!callbackResult(f, partition))
            return result_break;
        partition.pop_back();

        T lb = partition.empty() ? T(2) : partition.back();
        for (T i = isqrt(n); i >= lb; i--)
        {
            if (n % i == 0)
            {
                partition.push_back(i);
                if (!_enumerate(partition, n / i, f))
                    return result_break;
                partition.pop_back();
            }
        }
        return result_continue;
    }
};
} // namespace it
} // namespace euler
