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
/// Quadratic equation, with methods to solve mod primes and powers of primes. Precondition: gcd(a, b, c) = 1.
struct Quadratic
{
    i64 a, b, c;

    /// Solves ax^2 + bx + c == 0 (mod p). Precondition: gcd(a, b, c) = 1.
    [[nodiscard]] constexpr std::vector<i64> solveModPrime(i64 p) const
    {
        auto const a = mod(this->a, p);
        auto const b = mod(this->b, p);
        auto const c = mod(this->c, p);
        // If a == 0, this is actually linear.
        if (a == 0)
        {
            if (b == 0)
                return {};
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

    /// Solves ax^2 + bx + c == 0 (mod p^e). Precondition: gcd(a, b, c) = 1.
    [[nodiscard]] constexpr std::vector<i64> solveModPrimePower(i64 p, int e) const
    {
        if (p == 2)
        {
            // Brute force (TODO: do this properly at some point).
            i64 const q = 1_i64 << e;
            return it::range(0, q - 1).filter([&](i64 x) { return (a * x * x + b * x + c) % q == 0; }).to();
        }
        std::vector<std::pair<i64, int>> v = mapv(solveModPrime(p), [](i64 x) { return std::pair{x, 1}; });
        std::vector<i64> res;
        while (!v.empty())
        {
            auto const [x, i] = v.back();
            v.pop_back();
            if (i == e)
                res.push_back(x);
            else
            {
                // Hensel lifting.
                auto const q = pow(p, i);
                auto const num = mod(a * x * x + b * x + c, p * q);
                auto const den = 2 * a * x + b;
                if (den % p != 0)
                    v.emplace_back(mod(x - num * modInverse(den, p * q), p * q), i + 1);
                else if (num == 0)
                    for (i64 k = 0; k < p; ++k)
                        v.emplace_back(x + k * q, i + 1);
            }
        }
        std::ranges::sort(res);
        return res;
    }

    /// Solves ax^2 + bx + c == 0 (mod n). Precondition: gcd(a, b, c) = 1.
    [[nodiscard]] constexpr std::vector<i64> solveMod(i64 n) const
    {
        std::vector<i64> res{0};
        std::vector<i64> copy;
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
            swap(copy, res);
            res.clear();
            for (auto x : copy)
                for (auto y : sols)
                    res.push_back(crt(x, y, m, q));
            m *= q;
            return it::result_continue;
        });
        std::ranges::sort(res);
        return res;
    }
};
}; // namespace euler
