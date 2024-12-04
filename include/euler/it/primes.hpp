#pragma once

#include "base.hpp"
#include <primesieve/iterator.hpp>

inline namespace euler
{
namespace it
{
/// Enumerates primes, using the primesieve library.
class primes : public it_base
{
  public:
    using value_type = int64_t;

    constexpr primes(int64_t start = 2, int64_t stop = std::numeric_limits<int64_t>::max()) : _start(start), _stop(stop)
    {
    }

    template <std::invocable<value_type> Fun> constexpr result_t operator()(Fun f) const
    {
        if (_start > _stop)
            return result_continue;
        primesieve::iterator it(
            _start, _stop == std::numeric_limits<int64_t>::max() ? std::numeric_limits<uint64_t>::max() : _stop);
        for (int64_t p = it.next_prime(); p <= _stop; p = it.next_prime())
            if (!callbackResult(f, p))
                return result_break;
        return result_continue;
    }

  private:
    int64_t _start;
    int64_t _stop;
};
} // namespace it
} // namespace euler
