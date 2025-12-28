#pragma once

#include "algorithm.hpp"

namespace euler
{ /// Prime factorization class.
template <typename T = int64_t> class PF
{
    std::vector<PrimePower<T>> data_;

  public:
    using value_type = PrimePower<T>;
    using container_type = std::vector<value_type>;
    using reference = container_type::reference;
    using const_reference = container_type::const_reference;
    using iterator = container_type::iterator;
    using const_iterator = container_type::const_iterator;

    PF() = default;

    /// Initialize with a vector of pairs.
    /// Precondition: the vector must be in canonical form (first elements are primes in increasing order).
    constexpr PF(container_type pairs) : data_(std::move(pairs)) {}

    /// Initialize with a pair.
    constexpr PF(value_type pair) : data_({pair}) {}

    template <std::ranges::range Range>
    constexpr PF(std::from_range_t f, Range &&r) : data_(std::move(f), std::forward<Range>(r))
    {
    }

    /// Mutable accessor to data. Make sure invariants are preserved.
    [[nodiscard]] constexpr container_type &data() noexcept { return data_; }
    [[nodiscard]] constexpr const container_type &data() const noexcept { return data_; }
    [[nodiscard]] constexpr size_t size() const noexcept { return data_.size(); }
    [[nodiscard]] constexpr bool empty() const noexcept { return data_.empty(); }
    /// Mutable begin iterator. Make sure invariants are preserved.
    [[nodiscard]] constexpr auto begin() noexcept { return data_.begin(); }
    [[nodiscard]] constexpr auto begin() const noexcept { return data_.begin(); }
    /// Mutable end iterator. Make sure invariants are preserved.
    [[nodiscard]] constexpr auto end() noexcept { return data_.end(); }
    [[nodiscard]] constexpr auto end() const noexcept { return data_.end(); }
    [[nodiscard]] constexpr reference front() noexcept { return data_.front(); }
    [[nodiscard]] constexpr const_reference front() const noexcept { return data_.front(); }
    [[nodiscard]] constexpr reference back() noexcept { return data_.back(); }
    [[nodiscard]] constexpr const_reference back() const noexcept { return data_.back(); }

    [[nodiscard]] constexpr reference operator[](size_t i) noexcept { return data_[i]; }
    [[nodiscard]] constexpr const_reference operator[](size_t i) const noexcept { return data_[i]; }

    constexpr void pop_back() noexcept { data_.pop_back(); }
    constexpr void swap(PF &other) noexcept { data_.swap(other.data_); }
    constexpr void clear() noexcept { data_.clear(); }
    /// Dangerous! Make sure you know what you are doing.
    template <typename... Args> constexpr reference emplace_back(Args &&...args)
    {
        return data_.emplace_back(std::forward<Args>(args)...);
    }
    /// Dangerous! Make sure you know what you are doing.
    constexpr void push_back(value_type x) { data_.push_back(x); }
    /// Dangerous! Make sure you know what you are doing.
    constexpr iterator insert(const_iterator position, value_type x) { return data_.insert(position, x); }

    /// Returns the exponent of `p` in the prime factorization.
    constexpr int exponent(T p) const
    {
        if (auto const it = std::ranges::lower_bound(data_, p, std::ranges::less{}, [](auto &&t) { return t.first; });
            it != end() && it->first == p)
            return it->second;
        return 0;
    }

    /// Combines two prime factorizations using a binary operation `op` on the exponents.
    template <typename U, std::invocable<int, int> BinaryOp>
    constexpr PF &combineInPlace(const PF<U> &other, BinaryOp op)
    {
        auto it1 = data_.begin();
        auto it2 = other.begin();
        while (it1 != end() && it2 != other.end())
        {
            auto &[p1, e1] = *it1;
            auto [p2, e2] = *it2;
            if (p1 == p2)
            {
                e1 = op(e1, e2);
                if (e1 == 0)
                    it1 = data_.erase(it1);
                else
                    ++it1;
                ++it2;
            }
            else if (p1 < p2)
            {
                e1 = op(e1, 0);
                if (e1 == 0)
                    it1 = data_.erase(it1);
                else
                    ++it1;
            }
            else
            {
                auto e = op(0, e2);
                if (e != 0)
                    it1 = data_.emplace(it1, p2, e) + 1;
                ++it2;
            }
        }
        for (; it2 != other.end(); ++it2)
        {
            auto [p2, e2] = *it2;
            auto e = op(0, e2);
            if (e != 0)
                data_.emplace_back(p2, e);
        }
        return *this;
    }

    constexpr PF &powInPlace(int n) noexcept
    {
        for (auto &[p, e] : data_)
            e *= n;
        return *this;
    }

    [[nodiscard]] constexpr PF pow(this PF self, int n) noexcept
    {
        self.powInPlace(n);
        return self;
    }

    /// Checks if `this` divides `other`.
    [[nodiscard]] constexpr bool divides(const PF &other) const
    {
        auto it1 = begin();
        auto it2 = other.begin();
        while (it1 != end() && it2 != other.end())
        {
            if (it1->first == it2->first)
            {
                if (it1->second > it2->second)
                    return false;
                ++it1;
                ++it2;
            }
            else if (it1->first > it2->first)
                ++it2;
            else
                return false;
        }
        return it1 == end();
    }

    /// Evaluates the prime factorization to a number.
    template <typename Z = T> constexpr Z value() const
    {
        return product(data_, [](auto &&pe) { return euler::pow(Z(pe.first), pe.second); });
    }

    /// Counts the number of divisors from the factorization.
    template <typename Z = T> [[nodiscard]] constexpr Z countDivisors() const
    {
        return product(data_, [](auto &&pe) { return Z(pe.second + 1); });
    }

    /// Sums the divisors from the factorization.
    template <typename Z = T> [[nodiscard]] constexpr Z sumDivisors() const
    {
        return product(data_,
                       [](auto &&pe) { return Z(euler::pow(Z(pe.first), pe.second + 1) - 1) / Z(pe.first - 1); });
    }

    /// Computes Euler's totient function.
    template <typename Z = T> [[nodiscard]] constexpr Z totient() const
    {
        Z n = value<Z>();
        for (auto &&[p, _] : data_)
            n = n / p * (p - 1);
        return n;
    }

    /// Computes the radical.
    template <typename Z = T> [[nodiscard]] constexpr Z radical() const
    {
        Z rad(1);
        for (auto &&[p, _] : data_)
            rad *= p;
        return rad;
    }

    /// Returns whether the number is prime (has only one factor, with exponent 1).
    [[nodiscard]] constexpr bool isPrime() const { return data_.size() == 1 && data_[0].second == 1; }

    /// Multiplication.
    template <typename U> constexpr PF &operator*=(const PF<U> &other) { return combineInPlace(other, std::plus{}); }

    template <typename U> [[nodiscard]] constexpr friend PF operator*(PF left, const PF<U> &right)
    {
        left *= right;
        return left;
    }

    /// Division.
    template <typename U> constexpr PF &operator/=(const PF<U> &other) { return combineInPlace(other, std::minus{}); }

    template <typename U> [[nodiscard]] constexpr friend PF operator/(PF left, const PF<U> &right)
    {
        left /= right;
        return left;
    }

    /// LCM.
    template <typename U> constexpr PF &operator|=(const PF<U> &other) { return combineInPlace(other, maximum{}); }

    template <typename U> [[nodiscard]] constexpr friend PF operator|(PF left, const PF<U> &right)
    {
        left |= right;
        return left;
    }

    /// GCD.
    template <typename U> constexpr PF &operator&=(const PF<U> &other) { return combineInPlace(other, minimum{}); }

    template <typename U> [[nodiscard]] constexpr friend PF operator&(PF left, const PF<U> &right)
    {
        left &= right;
        return left;
    }

    /// Spaceship (3-way comparison) operator
    std::strong_ordering operator<=>(const PF &other) const = default;

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const PF &pf)
    {
        if (pf.empty())
            return o << '1';
        for (auto i = pf.begin(); i != pf.end(); ++i)
        {
            if (i != pf.begin())
                o << " * ";
            // Add parentheses if our type isn't integral, just to be safe
            if constexpr (integral2<T>)
                o << i->first;
            else
                o << '(' << i->first << ')';
            if (i->second != 1)
                o << '^' << i->second;
        }
        return o;
    }
};
} // namespace euler
