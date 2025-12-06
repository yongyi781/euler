#pragma once

#include "algorithm.hpp"
#include "decls.hpp"
#include "it/base.hpp"
#include "it/factor.hpp"
#include "literals.hpp"
#include "math.hpp"
#include "modular_arithmetic.hpp"
#include "types.hpp"

namespace euler
{
/// Quadratic equation, with methods to solve mod primes and powers of primes.
struct Quadratic
{
    i64 a, b, c;

    /// Solves ax^2 + bx + c == 0 (mod p).
    [[nodiscard]] constexpr std::vector<i64> solveModPrime(i64 p) const
    {
        if (p == 2)
        {
            std::vector<i64> res;
            if (c % 2 == 0)
                res.push_back(0);
            if ((a + b + c) % 2 == 0)
                res.push_back(1);
            return res;
        }
        auto const a = mod(this->a, p);
        auto const b = mod(this->b, p);
        auto const c = mod(this->c, p);
        // If a == 0, this is actually linear.
        if (a == 0)
        {
            if (b == 0)
                return c == 0 ? range(0, p - 1) : std::vector<i64>{};
            return {mod(-c * modInverse(b, p), p)};
        }
        auto const d = mod(b * b - 4 * a * c, p);
        auto const inv_2a = modInverse(2 * a, p);
        if (d == 0)
            return {mod(-b * inv_2a, p)};
        auto const res = sqrtModp(d, p);
        if (res == -1)
            return {};
        auto const x1 = mod((-b + res) * inv_2a, p), x2 = mod((-b - res) * inv_2a, p);
        if (x1 < x2)
            return {x1, x2};
        return {x2, x1};
    }

    /// Solves ax^2 + bx + c == 0 (mod p^e).
    [[nodiscard]] constexpr std::vector<i64> solveModPrimePower(i64 p, int e) const
    {
        auto sols = solveModPrime(p);
        if (sols.empty())
            return {};
        std::vector<i64> next_sols;
        i64 q = p;
        for (; --e > 0; q *= p)
        {
            // Hensel lifting.
            next_sols.clear();
            auto const q2 = q * p;
            for (auto &&x : sols)
            {
                auto const num = mod(a * x % q2 * x % q2 + b * x + c, q2);
                auto const den = 2 * a * x + b;
                if (den % p != 0)
                    next_sols.push_back(mod(x - num * modInverse(den, q2), q2));
                else if (num == 0)
                    for (i64 k = 0; k < p; ++k)
                        next_sols.push_back(x + k * q);
            }
            if (next_sols.empty())
                return {};
            swap(sols, next_sols);
        }
        std::ranges::sort(sols);
        return sols;
    }

    /// Solves ax^2 + bx + c == 0 (mod n).
    [[nodiscard]] constexpr std::vector<i64> solveMod(i64 n) const
    {
        std::vector<i64> res{0};
        std::vector<i64> next_res;
        i64 m = 1;
        it::factor{n}([&](auto &&pe) {
            auto &&[p, e] = pe;
            i64 const q = pow(p, e);
            auto const sols = solveModPrimePower(p, e);
            if (sols.empty())
            {
                res.clear();
                return it::result_break;
            }
            next_res.clear();
            for (auto x : res)
                for (auto y : sols)
                    next_res.push_back(crt(x, y, m, q));
            m *= q;
            swap(res, next_res);
            return it::result_continue;
        });
        std::ranges::sort(res);
        return res;
    }
};
}; // namespace euler
