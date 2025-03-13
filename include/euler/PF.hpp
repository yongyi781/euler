#pragma once

#include "algorithm.hpp"

inline namespace euler
{ /// Prime factorization class.
template <typename T = int64_t> class PF
{
  public:
    using value_type = PrimePower<T>;
    using container_type = std::vector<value_type>;

    PF() = default;

    /// Initialize with a vector of pairs.
    /// Precondition: the vector must be in canonical form (first elements are primes in increasing order).
    constexpr PF(container_type pairs) : _data(std::move(pairs)) {}

    /// Initialize with a pair.
    constexpr PF(value_type pair) : _data({pair}) {}

    /// Mutable accessor to data. Make sure invariants are preserved.
    [[nodiscard]] constexpr container_type &data() noexcept { return _data; }
    [[nodiscard]] constexpr const container_type &data() const noexcept { return _data; }
    [[nodiscard]] constexpr size_t size() const noexcept { return _data.size(); }
    [[nodiscard]] constexpr bool empty() const noexcept { return _data.empty(); }
    /// Mutable begin iterator. Make sure invariants are preserved.
    [[nodiscard]] constexpr auto begin() noexcept { return _data.begin(); }
    [[nodiscard]] constexpr auto begin() const noexcept { return _data.begin(); }
    /// Mutable end iterator. Make sure invariants are preserved.
    [[nodiscard]] constexpr auto end() noexcept { return _data.end(); }
    [[nodiscard]] constexpr auto end() const noexcept { return _data.end(); }
    [[nodiscard]] constexpr value_type &front() noexcept { return _data.front(); }
    [[nodiscard]] constexpr const value_type &front() const noexcept { return _data.front(); }
    [[nodiscard]] constexpr value_type &back() noexcept { return _data.back(); }
    [[nodiscard]] constexpr const value_type &back() const noexcept { return _data.back(); }

    [[nodiscard]] constexpr value_type &operator[](size_t i) noexcept { return _data[i]; }
    [[nodiscard]] constexpr const value_type &operator[](size_t i) const noexcept { return _data[i]; }

    constexpr void pop_back() noexcept { _data.pop_back(); }
    constexpr void swap(PF &other) noexcept { _data.swap(other._data); }
    constexpr void clear() noexcept { _data.clear(); }

    /// Combines two prime factorizations using a binary operation `op` on the exponents.
    template <typename U, std::invocable<int, int> BinaryOp>
    constexpr PF &combineInPlace(const PF<U> &other, BinaryOp op)
    {
        auto it1 = _data.begin();
        auto it2 = other.begin();
        while (it1 != end() && it2 != other.end())
        {
            auto &[p1, e1] = *it1;
            auto [p2, e2] = *it2;
            if (p1 == p2)
            {
                e1 = op(e1, e2);
                if (e1 == 0)
                    it1 = _data.erase(it1);
                else
                    ++it1;
                ++it2;
            }
            else if (p1 < p2)
            {
                e1 = op(e1, 0);
                if (e1 == 0)
                    it1 = _data.erase(it1);
                else
                    ++it1;
            }
            else
            {
                auto e = op(0, e2);
                if (e != 0)
                    it1 = _data.emplace(it1, p2, e) + 1;
                ++it2;
            }
        }
        for (; it2 != other.end(); ++it2)
        {
            auto [p2, e2] = *it2;
            auto e = op(0, e2);
            if (e != 0)
                _data.emplace_back(p2, e);
        }
        return *this;
    }

    constexpr PF &powInPlace(int n) noexcept
    {
        for (auto &[p, e] : _data)
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
        return product(_data, [](auto &&pe) { return euler::pow(Z(pe.first), pe.second); });
    }

    /// Counts the number of divisors from the factorization.
    template <typename Z = T> [[nodiscard]] constexpr Z countDivisors() const
    {
        return product(_data, [](auto &&pe) { return Z(pe.second + 1); });
    }

    /// Sums the divisors from the factorization.
    template <typename Z = T> [[nodiscard]] constexpr Z sumDivisors() const
    {
        return product(_data,
                       [](auto &&pe) { return Z(euler::pow(Z(pe.first), pe.second + 1) - 1) / Z(pe.first - 1); });
    }

    /// Computes Euler's totient function.
    template <typename Z = T> [[nodiscard]] constexpr Z totient() const
    {
        Z n = value<Z>();
        for (auto &&[p, _] : _data)
            n = n / p * (p - 1);
        return n;
    }

    /// Computes the radical.
    template <typename Z = T> [[nodiscard]] constexpr Z radical() const
    {
        Z rad{1};
        for (auto &&[p, _] : _data)
            rad *= p;
        return rad;
    }

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
            o << i->first;
            if (i->second != 1)
                o << '^' << i->second;
        }
        return o;
    }

  private:
    container_type _data;
};
} // namespace euler
