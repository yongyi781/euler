#pragma once

#include <iterator>

namespace euler
{
// Removes the dependency on boost::counting_iterator.
template <typename T> struct counting_iterator
{
    using iterator_category = std::random_access_iterator_tag;
    using iterator_concept = std::random_access_iterator_tag;
    using difference_type = std::iter_difference_t<T>;
    using pointer = const T *;
    using reference = const T &;
    using value_type = T;

    T value;

    counting_iterator() = default;
    constexpr counting_iterator(T value) : value(std::move(value)) {}

    constexpr T operator*() const { return value; }

    constexpr counting_iterator &operator++()
    {
        ++value;
        return *this;
    }

    constexpr counting_iterator operator++(int)
    {
        counting_iterator tmp = *this;
        ++(*this);
        return tmp;
    }

    constexpr counting_iterator &operator+=(T i)
    {
        value += i;
        return *this;
    }

    constexpr counting_iterator operator+(const difference_type other) const { return {T(value + other)}; }
    friend constexpr counting_iterator operator+(const difference_type value, const counting_iterator &other)
    {
        return {T(other + value)};
    }

    constexpr counting_iterator &operator--()
    {
        --value;
        return *this;
    }
    constexpr counting_iterator operator--(int)
    {
        counting_iterator tmp = *this;
        --(*this);
        return tmp;
    }
    constexpr counting_iterator &operator-=(T i)
    {
        value -= i;
        return *this;
    }
    constexpr difference_type operator-(const counting_iterator &other) const { return value - other.value; }

    constexpr counting_iterator operator-(const difference_type other) const { return value - other; }
    friend constexpr counting_iterator operator-(const difference_type value, const counting_iterator &other)
    {
        return other - value;
    }

    constexpr T operator[](difference_type i) const { return value + i; }

    constexpr std::strong_ordering operator<=>(const counting_iterator &other) const
    {
        if (value == other.value)
            return std::strong_ordering::equal;
        return value < other.value ? std::strong_ordering::less : std::strong_ordering::greater;
    }

    constexpr bool operator==(const counting_iterator &other) const { return value == other.value; }
    constexpr bool operator!=(const counting_iterator &other) const { return value != other.value; }
};
} // namespace euler
