#pragma once

#include "Quadratic.hpp"
#include "cfrac.hpp"
#include "literals.hpp"
#include "matrix.hpp"
#include "types.hpp"

inline namespace euler
{
/// Quadratic diophantine equation.
struct QD
{
    /// The continued fraction for `(a + √b) / c`, for one period only.
    template <integral2 T> class cfrac_quadratic_one_period : public cfrac::cfrac_base
    {
      public:
        using value_type = T;

        T a, b, c;

        cfrac_quadratic_one_period() = default;
        constexpr cfrac_quadratic_one_period(T a, T b, T c) : a(std::move(a)), b(std::move(b)), c(std::move(c))
        {
            assert(this->b >= 0);
        }

        /// Enumerates the terms in the square root continued fraction.
        template <typename Fun> it::result_t operator()(Fun f) const
        {
            using DT = double_integer_t<T>;
            // Idea: each remainder is of the form (x + c√b) / d, starting with x = a*c and d = c².
            T const fl = isqrt(c * c * b) * (c > 0 ? 1 : -1);
            DT x = a * c;
            T d = c * c;
            // k = floor.
            T k = floorDiv(a * c + fl, c * c);
            if (!it::callbackResult(f, k))
                return it::result_break;
            boost::unordered_flat_set<std::pair<DT, T>> seen{{x, d}};
            while (true)
            {
                x = k * d - x;
                d = T((c * c * b - x * x) / d);
                if (d == 0 || !seen.emplace(x, d).second)
                    break;
                k = T(floorDiv(x + fl, d));
                if (!it::callbackResult(f, k))
                    return it::result_break;
            }
            return it::result_continue;
        }
    };

    /// Represents the affine linear transformation `(x, y) -> (p * x + q * y + k, r * x + s * y + l)`.
    struct affine_transform
    {
        i64 p, q, r, s, k, l;

        template <typename T = i64> [[nodiscard]] constexpr Vector<T, 2> operator()(T x, T y) const noexcept
        {
            return Vector<T, 2>{p * x + q * y + k, r * x + s * y + l};
        }

        template <typename T = i64> [[nodiscard]] constexpr Vector<T, 2> operator()(Vector<T, 2> v) const noexcept
        {
            return (*this)(v[0], v[1]);
        }

        /// Precondition: ps - qr = 1.
        [[nodiscard]] constexpr affine_transform inverse() const noexcept
        {
            return {.p = s, .q = -q, .r = -r, .s = p, .k = q * l - s * k, .l = r * k - p * l};
        }

        template <typename CharT, typename Traits>
        friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o,
                                                             const affine_transform &a)
        {
            return o << '(' << a.p << "x + " << a.q << "y + " << a.k << ", " << a.r << "x + " << a.s << "y + " << a.l
                     << ')';
        }
    };

    i64 a, b, c, d, e, f;

    QD() = default;
    constexpr QD(i64 a, i64 b, i64 c, i64 d, i64 e, i64 f) : a(a), b(b), c(c), d(d), e(e), f(f)
    {
        i64 const g = gcd({this->a, this->b, this->c, this->d, this->e, this->f}) * (f > 0 ? -1 : 1);
        this->a /= g;
        this->b /= g;
        this->c /= g;
        this->d /= g;
        this->e /= g;
        this->f /= g;
    }

    /// Returns the discriminant.
    [[nodiscard]] constexpr i64 disc() const noexcept { return b * b - 4 * a * c; }

    /// Evaluates the quadratic function at (x, y).
    template <integral2 T> [[nodiscard]] constexpr T operator()(const T &x, const T &y)
    {
        return a * x * x + b * x * y + c * y * y + d * x + e * y + f;
    }

    /// Fundamental solutions to `ax^2 + bxy + cy^2 + dx + ey + f == 0`.
    template <integral2 T = i64> [[nodiscard]] std::vector<Vector<T, 2>> solve() const
    {
        i64 const D = disc();
        assert(D > 0 && !isSquare(D) && "Non-hyperbolic QDs not implemented");

        if (d == 0 && e == 0 && std::gcd(a, f) == 1)
            return solve_af_coprime<T>(0, 0, 1, std::array<i64, 4>{1, 0, 0, 1});

        auto const alpha = 2 * c * d - b * e;
        auto const beta = 2 * a * e - b * d;
        i64 const n = -D * (a * e * e - b * e * d + c * d * d + f * D);
        if (n == 0)
            return {Vector<T, 2>{}};

        QD q{a, b, c, 0, 0, -n};
        auto const aa = q.a, bb = q.b, cc = q.c;
        std::array<i64, 4> coeffs{1, 0, 0, 1};
        if (std::gcd(aa, n) != 1)
        {
            // Must apply a unimodular transformation.
            i64 m = 0;
            while (true)
            {
                ++m;
                q.a = aa * m * m + bb * m + cc;
                if (std::gcd(q.a, n) == 1)
                {
                    q.b = 2 * (aa * m * m + (bb - aa) * m + cc) - bb;
                    q.c = aa * m * m + (bb - 2 * aa) * m + aa - bb + cc;
                    coeffs = {m, m - 1, 1, 1};
                    break;
                }

                q.a = aa + bb * (m - 1) + cc * (m - 1) * (m - 1);
                if (std::gcd(q.a, n) == 1)
                {
                    q.b = 2 * (cc * m * m + (bb - cc) * m + aa) - bb;
                    q.c = aa + bb * m + cc * m * m;
                    coeffs = {1, 1, m - 1, m};
                    break;
                }
            }
        }
        return q.solve_af_coprime<T>(alpha, beta, D, coeffs);
    }

    /// Returns all integral solutions to `ax^2 + bxy + cy^2 + dx + ey + f == 0` within the specified range.
    template <integral2 T> [[nodiscard]] std::vector<Vector<T, 2>> solve(T xmin, T xmax, T ymin, T ymax)
    {
        using std::abs;

        auto sols = solve();
        auto const rec = recurrence();
        if (sols.size() == 1 && sols[0][0] == 0 && sols[0][1] == 0 && rec.k == 0 && rec.l == 0)
            return sols;
        auto const inv = rec.inverse();
        T const xlimit = std::max({abs(xmin), abs(xmax)});
        T const ylimit = std::max({abs(ymin), abs(ymax)});
        std::vector<Vector<T, 2>> res;
        for (auto &&v0 : sols)
        {
            Vector<T, 2> v{v0[0], v0[1]};
            while (abs(v[0]) <= xlimit && abs(v[1]) <= ylimit)
            {
                if (v[0] >= xmin && v[0] <= xmax && v[1] >= ymin && v[1] <= ymax)
                    res.push_back(v);
                v = rec(v);
            }
            v = {v0[0], v0[1]};
            v = inv(v);
            while (abs(v[0]) <= xlimit && abs(v[1]) <= ylimit)
            {
                if (v[0] >= xmin && v[0] <= xmax && v[1] >= ymin && v[1] <= ymax)
                    res.push_back(v);
                v = inv(v);
            }
        }

        std::ranges::sort(res);
        auto const u = std::ranges::unique(res);
        res.erase(u.begin(), u.end());
        return res;
    }

    /// Gets the forward recurrence for solutions to this quadratic diophantine equation.
    template <integral2 T = i64> [[nodiscard]] affine_transform recurrence() const
    {
        i64 const D = disc();
        auto const res = cfrac::quadratic(-b, b * b - 4 * a * c, 2_i64)
                             .convergents()
                             .filterMap([&](auto &&t) -> std::optional<affine_transform> {
                                 auto &&[x, y] = t;
                                 if (x * x + b * x * y + a * c * y * y != 1)
                                     return std::nullopt;
                                 for (auto &&[r, s] : {std::pair{x, y}, std::pair{-x, -y}})
                                 {
                                     i64 const P = r;
                                     i64 const Q = -c * s;
                                     i64 const R = a * s;
                                     i64 const S = r + b * s;
                                     i64 K = -c * d * (P + S - 2) - e * (b - b * r - 2 * a * c * s);
                                     i64 L = -d * (b - b * r - 2 * a * c * s) - a * e * (P + S - 2);
                                     if (K % D == 0 && L % D == 0)
                                         return affine_transform{P, Q, R, S, K / D, L / D + d * s};
                                 }
                                 return std::nullopt;
                             })
                             .first();
        return *res;
    }

    template <typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &o, const QD &q)
    {
        return o << "QD: " << std::vector{q.a, q.b, q.c, q.d, q.e, q.f} << " (disc = " << q.disc() << ")";
    }

  private:
    /// Primitive undamental solutions to `ax^2 + bxy + cy^2 + f == 0`, subject to a predicate.
    /// Precondition: `gcd(a, b, c) = 1` and `gcd(a, f) = 1`.
    template <integral2 T = i64, typename Pred>
    [[nodiscard]] std::vector<Vector<T, 2>> solvePrimitive_af_coprime(Pred pred) const
    {
        std::vector<Vector<T, 2>> res;
        i64 const n = -f;
        auto const proc = [&](i64 t, i64 p, i64 q, i64 r, const T &k, const T &y) -> bool {
            bool found = false;
            if (p * y * y + q * y * k + r * k * k == 1)
            {
                T const x = t * y - n * k;
                if (pred(x, y))
                {
                    found = true;
                    res.emplace_back(x, y);
                }
                if (pred(-x, -y))
                {
                    found = true;
                    res.emplace_back(-x, -y);
                }
            }
            return !found;
        };
        for (i64 t : Quadratic(a, b, c).solveMod(n))
        {
            if (t > n / 2)
                t -= n;
            i64 const p = (a * t * t + b * t + c) / n;
            if (p == 1)
            {
                bool found = false;
                if (pred(t, 1))
                {
                    found = true;
                    res.emplace_back(t, 1);
                }
                if (pred(-t, -1))
                {
                    found = true;
                    res.emplace_back(-t, -1);
                }
                if (found)
                    continue;
            }
            i64 const q = -(2 * a * t + b);
            i64 const r = a * n;
            if (gcd({p, q, r}) != 1)
                continue;

            cfrac_quadratic_one_period<T>(-q, disc(), 2 * r).template convergents<T>()([&](auto &&ky) {
                return proc(t, p, q, r, ky.first, ky.second);
            });
            cfrac_quadratic_one_period<T>(q, disc(), -2 * r).template convergents<T>()([&](auto &&ky) {
                return proc(t, p, q, r, ky.first, ky.second);
            });
        }
        return res;
    }

    /// Fundamental solutions to ax^2 + bxy + cy^2 + f == 0 and gcd(x, y) = 1.
    /// Precondition: gcd(a, b, c) = 1 and gcd(a, f) = 1.
    /// Here, D is not the discriminant but rather the denominator of the Legendre transformation.
    template <integral2 T = i64>
    [[nodiscard]] std::vector<Vector<T, 2>> solve_af_coprime(i64 alpha, i64 beta, i64 D,
                                                             std::array<i64, 4> coeffs) const
    {
        std::vector<Vector<T, 2>> res;
        // Group solutions by g = gcd(x, y).
        QD q = *this;
        for (auto const g : divisors(sqfreeDecompose(-f).first))
        {
            q.f = f / (g * g);
            for (auto const [u, v] : q.solvePrimitive_af_coprime<T>([&](auto &&u, auto &&v) {
                     return (g * (coeffs[0] * u + coeffs[1] * v) + alpha) % D == 0 &&
                            (g * (coeffs[2] * u + coeffs[3] * v) + beta) % D == 0;
                 }))
            {
                res.emplace_back((g * (coeffs[0] * u + coeffs[1] * v) + alpha) / D,
                                 (g * (coeffs[2] * u + coeffs[3] * v) + beta) / D);
            }
        }
        std::ranges::sort(res);
        auto const u = std::ranges::unique(res);
        res.erase(u.begin(), u.end());
        return res;
    }
};
} // namespace euler
