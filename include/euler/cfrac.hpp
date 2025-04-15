#pragma once

#include "decls.hpp"
#include "io.hpp"
#include "it.hpp"

inline namespace euler
{
namespace cfrac
{
template <typename T, it::enumerable CFrac> class convergents_t;

/// The base class for continued fractions. Don't use this class directly.
class cfrac_base : public it::it_base
{
  public:
    /// Enumerates the convergents.
    template <typename T = int64_t, typename Self>
    constexpr convergents_t<T, std::decay_t<Self>> convergents(this Self &&self)
    {
        return {std::forward<Self>(self)};
    }
};

/// Makes a continued fraction out of a list of numbers.
template <std::ranges::view V> class from_terms : public cfrac_base
{
  public:
    using value_type = std::ranges::range_value_t<V>;

    from_terms() = default;
    constexpr explicit from_terms(V base) : _terms(std::move(base)) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        for (auto &&x : _terms)
            if (!it::callbackResult(f, x))
                return it::result_break;
        return it::result_continue;
    }

  private:
    V _terms;
};

template <std::ranges::range Range> from_terms(Range &&) -> from_terms<std::views::all_t<Range>>;

/// The continued fraction for the square root of a nonnegative integer. Use this instead of
/// sqrtAsPeriodic() if you call this in a loop, since this one does not allocate.
template <integral2 T> class sqrt : public cfrac_base
{
  public:
    using value_type = T;

    T radicand;

    sqrt() = default;
    constexpr sqrt(T radicand) : radicand(std::move(radicand)) { assert(this->radicand >= 0); }

    /// Enumerates the terms in the square root continued fraction.
    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        T const fl = isqrt(radicand);
        T x = 0;
        T d = 1;
        // k = floor.
        T k = fl;
        if (!it::callbackResult(f, k))
            return it::result_break;
        while (true)
        {
            x = d * k - x;
            d = (radicand - x * x) / d;
            if (d == 0)
                break;
            k = (fl + x) / d;
            if (!it::callbackResult(f, k))
                return it::result_break;
        }
        return it::result_continue;
    }
};

/// The continued fraction for `(a + √b) / c`.
template <integral2 T> class quadratic : public cfrac_base
{
  public:
    using value_type = T;

    T a, b, c;

    quadratic() = default;
    constexpr quadratic(T a, T b, T c) : a(std::move(a)), b(std::move(b)), c(std::move(c)) { assert(this->b >= 0); }

    /// Enumerates the terms in the square root continued fraction.
    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        // Idea: each remainder is of the form (x + c√b) / d, starting with x = a*c and d = c².
        T const fl = isqrt(c * c * b);
        T x = a * c;
        T d = c * c;
        // k = floor.
        T k = (a * c + fl) / (c * c);
        if (!it::callbackResult(f, k))
            return it::result_break;
        while (true)
        {
            x = k * d - x;
            d = (c * c * b - x * x) / d;
            if (d == 0)
                break;
            k = (x + fl) / d;
            if (!it::callbackResult(f, k))
                return it::result_break;
        }
        return it::result_continue;
    }
};

/// A periodic continued fraction.
template <integral2 T> class periodic : public cfrac_base
{
  public:
    using value_type = T;

    periodic() = default;
    constexpr periodic(std::vector<T> leadingTerms, std::vector<T> periodicTerms)
        : _leadingTerms(std::move(leadingTerms)), _periodicTerms(std::move(periodicTerms))
    {
    }
    [[nodiscard]] constexpr const std::vector<T> &leadingTerms() const { return _leadingTerms; }
    [[nodiscard]] constexpr const std::vector<T> &periodicTerms() const { return _periodicTerms; }

    /// Enumerates the terms in the periodic continued fraction.
    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        for (auto &&x : _leadingTerms)
            if (!it::callbackResult(f, x))
                return it::result_break;
        if (!_periodicTerms.empty())
            for (size_t i = 0;; i = (i + 1) % _periodicTerms.size())
                if (!it::callbackResult(f, _periodicTerms[i % _periodicTerms.size()]))
                    return it::result_break;
        return it::result_continue;
    }

    [[nodiscard]] constexpr std::string str() const
    {
        std::stringstream ss;
        ss << "[";
        for (auto &&x : _leadingTerms)
            ss << x << ", ";
        io::print(_periodicTerms, std::numeric_limits<size_t>::max(), ss);
        ss << "]";
        return std::move(ss).str();
    }

  private:
    std::vector<T> _leadingTerms;
    std::vector<T> _periodicTerms;
};

/// The continued fraction of Euler's constant e.
struct e_t : public cfrac_base
{
    using value_type = int64_t;

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        if (!it::callbackResult(f, (int64_t)2))
            return it::result_break;
        for (int64_t n = 2;; n += 2)
        {
            if (!it::callbackResult(f, (int64_t)1))
                return it::result_break;
            if (!it::callbackResult(f, n))
                return it::result_break;
            if (!it::callbackResult(f, (int64_t)1))
                return it::result_break;
        }
    }
};

/// The continued fraction of Euler's constant e.
constexpr e_t e;

/// The continued fraction of a rational number.
template <integral2 T> class rational : public cfrac_base
{
  public:
    using value_type = T;

    rational() = default;
    constexpr rational(T numerator, T denominator) : numerator(numerator), denominator(denominator) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        T n = numerator;
        T d = denominator;
        while (d != 0)
        {
            T x = n / d;
            if (!it::callbackResult(f, x))
                return it::result_break;
            n -= d * x;
            std::tie(n, d) = std::tuple{d, n};
        }
        return it::result_continue;
    }

  private:
    T numerator, denominator;
};

/// The continued fraction of a floating-point number.
template <integral2 Z = int64_t, typename T = double> class floating_point : public cfrac_base
{
  public:
    using value_type = Z;

    floating_point() = default;
    constexpr explicit floating_point(T value) : _value(value) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        using std::floor;

        T value = _value;
        while (true)
        {
            Z x = (Z)floor(value);
            if (!it::callbackResult(f, x))
                return it::result_break;
            value -= x;
            if (value == 0)
                break;
            value = T(1) / value;
            if (std::numeric_limits<Z>::max() > 0 && value > (T)std::numeric_limits<Z>::max())
                break;
        }
        return it::result_continue;
    }

  private:
    T _value;
};

/// Enumerates the convergents of the given continued fraction.
template <typename T, it::enumerable CFrac> class convergents_t : public it::it_base
{
  public:
    using value_type = std::pair<T, T>;

    convergents_t() = default;
    constexpr convergents_t(CFrac cfrac) : _cfrac(std::move(cfrac)) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        T h = 1;
        T k = 0;
        T h0 = 0;
        T k0 = 1;
        return _cfrac([&](auto &&x) {
            T h2 = x * h + h0;
            T k2 = x * k + k0;
            std::tie(h, k, h0, k0) = std::tuple{h2, k2, h, k};
            return it::callbackResult(f, std::pair{h, k});
        });
    }

  private:
    CFrac _cfrac;
};

/// Makes a periodic continued fraction for sqrt(d). This one allocates but should be faster
/// when constructing convergents.
template <integral2 T> constexpr periodic<T> sqrtAsPeriodic(const T &d)
{
    assert(d >= 0);
    std::vector<T> periodicTerms;
    T a0 = isqrt(d);
    T a = a0;
    T b = 0;
    T c = 1;
    while (a != 2 * a0)
    {
        b = c * a - b;
        c = (d - b * b) / c;
        if (c == 0)
            break;
        a = (a0 + b) / c;
        periodicTerms.push_back(a);
    }
    return periodic<T>({a0}, std::move(periodicTerms));
}

/// Enumerates solutions to the Pell equation x^2 - d*y^2 = 1.
template <integral2 T> class pell_solutions : public it::it_base
{
  public:
    using value_type = std::pair<T, T>;

    pell_solutions() = default;
    constexpr explicit pell_solutions(T d) : _d(std::move(d)) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        if (isSquare(_d))
            return it::result_break;
        return sqrt(_d).template convergents<T>()([&](auto &&r) {
            auto &&[x, y] = r;
            if (x * x - _d * y * y != 1)
                return it::result_continue;
            return it::callbackResult(f, std::forward<decltype(r)>(r));
        });
    }

  private:
    T _d;
};
} // namespace cfrac

/// Finds the best rational approximation to a floating-point number with denominator less than `denomBound`.
template <typename R, integral2 Z> std::pair<Z, Z> nearbyRational(R value, Z denomBound)
{
    using std::abs;

    assert(denomBound > 0);

    Z h0 = 0;
    Z k0 = 1;
    Z h1 = 1;
    Z k1 = 0;
    cfrac::floating_point{value}([&](const Z &a) {
        Z h2 = a * h1 + h0;
        Z k2 = a * k1 + k0;
        if (k2 > denomBound)
        {
            Z a2 = (denomBound - k0) / k1;
            if (a2 >= (a + 1) / 2)
            {
                h2 = a2 * h1 + h0;
                k2 = a2 * k1 + k0;
                if (abs((R)h2 / k2 - value) < abs((R)h1 / k1 - value))
                    std::tie(h0, h1, k0, k1) = std::tie(h1, h2, k1, k2);
            }
            return false;
        }
        std::tie(h0, h1, k0, k1) = std::tie(h1, h2, k1, k2);
        return true;
    });
    return {h1, k1};
}
} // namespace euler
