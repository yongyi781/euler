#pragma once

#include <boost/container_hash/hash.hpp>

namespace boost
{
template <class T> class rational;

template <typename T> constexpr size_t hash_value(const boost::rational<T> &r)
{
    size_t seed = 0;
    boost::hash_combine(seed, r.numerator());
    boost::hash_combine(seed, r.denominator());
    return seed;
}
} // namespace boost
