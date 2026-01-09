#pragma once

#include "euler/math.hpp"
#include "prime.hpp"

namespace euler
{
/// A class to compute binomial coefficients mod prime powers.
template <u64 P, int E> class BinomialModPrimePower
{
    static constexpr u64 PE = pow(P, E);

    // p^k
    std::vector<u64> pPower;

    // First kind
    std::vector<std::vector<u64>> stirling;

    // (p-1)!^k
    std::vector<u64> pstirling;

    // factorials and its inverse, but without the factor p
    std::vector<u64> fact, ifact;
    // p-adic valuation of factorial
    std::vector<int> pfact;

    std::vector<u64> lagrangeCoeff;

    std::vector<u64> prods;

    // Computes (np+1)(np+2) ... (np+m) % p^e
    // with m < p
    [[nodiscard]] u64 risingFactorial(u64 n, u64 m) const
    {
        u64 ret = 0;
        u64 pn = 1;
        int const ep = (int)(E < P ? E : P);
        for (int j = 0; j < ep; ++j)
        {
            ret = (u64)(((__int128)stirling[m + 1][j + 1] * pn + ret) % PE);
            pn = (u64)((__int128)pn * P * n % PE);
        }

        return ret;
    }

    [[nodiscard]] u64 lagrangeInterpolate(u64 ndp) const
    {
        int const len = E * 2 - 1;
        u64 ptot = 0;

        std::vector pfactorsnum(len, 0);
        std::vector prenum(len, (u64)0);
        std::vector sufnum(len, (u64)0);

        for (int i = 0; i < len; ++i)
        {
            u64 num = ndp - i;
            pfactorsnum[i] = 0;
            while (num % P == 0)
            {
                num /= P;
                ++pfactorsnum[i];
            }
            ptot += pfactorsnum[i];
            prenum[i] = sufnum[i] = num;

            if (i > 0)
                prenum[i] = (u64)((__int128)prenum[i - 1] * prenum[i] % PE);
        }

        for (int i = len - 2; i >= 0; --i)
            sufnum[i] = (u64)((__int128)sufnum[i] * sufnum[i + 1] % PE);

        u64 sum = 0;
        for (int j = 0; j < len; ++j)
        {
            int const j2 = len - 1 - j;
            u64 const pfactor = ptot - pfactorsnum[j] - pfact[j] - pfact[j2];

            if (pfactor >= E)
                continue;

            auto numerator = (u64)((__int128)(j > 0 ? prenum[j - 1] : 1) * (j < len - 1 ? sufnum[j + 1] : 1) % PE);

            sum = (u64)(((__int128)numerator * lagrangeCoeff[j] % PE * pPower[pfactor] + sum) % PE);
        }

        return sum;
    }

  public:
    // O(p*min(e, p) + e^2)
    BinomialModPrimePower() : pPower(powers(P, E)), fact(E * 2 - 1), ifact(E * 2 - 1), pfact(E * 2 - 1)
    {
        int const ep = (int)(E < P ? E : P);
        stirling = std::vector<std::vector<u64>>(P + 1, std::vector<u64>(ep + 1, 0));
        stirling[0][0] = 1;
        for (int i = 1; i <= P; ++i)
            for (int j = 1; j <= ep; ++j)
                stirling[i][j] = (u64)(((u128)(i - 1) * stirling[i - 1][j] + stirling[i - 1][j - 1]) % PE);

        prods.resize(E * 2 - 1);
        u64 prod = 1;
        u64 const invStirling = modInverse(stirling[P][1], PE);

        for (int i = 0; i <= E * 2 - 3; ++i)
        {
            prods[i] = prod;
            prod = (u64)((u128)prod * risingFactorial(i, P - 1) % PE * invStirling % PE);
        }
        prods[E * 2 - 2] = prod;

        pstirling.resize(E);
        pstirling[0] = 1;
        for (int i = 1; i < E; ++i)
            pstirling[i] = (u64)((u128)pstirling[i - 1] * stirling[P][1] % PE);

        int const len = E * 2 - 1;
        fact[0] = 1;
        pfact[0] = 0;
        for (int i = 1; i < len; ++i)
        {
            u64 num = i;
            pfact[i] = 0;
            while (num % P == 0)
            {
                ++pfact[i];
                num /= P;
            }

            ifact[i - 1] = num;
            pfact[i] += pfact[i - 1];
            fact[i] = (u64)((__int128)fact[i - 1] * num % PE);
        }

        ifact[len - 1] = modInverse(fact[len - 1], PE);
        for (int i = len - 2; i >= 0; --i)
            ifact[i] = (u64)((__int128)ifact[i + 1] * ifact[i] % PE);

        lagrangeCoeff.resize(len);
        for (int i = 0; i < len; ++i)
        {
            auto denominator = (u64)((__int128)ifact[i] * ifact[len - 1 - i] % PE);
            if (((len - 1 - i) & 1) != 0)
                denominator = PE - denominator;

            lagrangeCoeff[i] = (u64)((__int128)denominator * prods[i] % PE);
        }
    }

    [[nodiscard]] constexpr u64 p() const { return P; }

    [[nodiscard]] constexpr u64 e() const { return E; }

    [[nodiscard]] constexpr u64 pe() const { return PE; }

    // n! % p^e, but p, 2p, ... is not multiplied
    [[nodiscard]] u64 factorialWithoutP(u64 n) const
    {
        u64 const ndp = n / P;
        u64 ret = risingFactorial(ndp, n % P);

        if (ndp <= E * 2 - 2)
            return (u64)((__int128)ret * prods[ndp] % PE);

        ret = (u64)((__int128)ret * lagrangeInterpolate(ndp) % PE);
        return ret;
    }

    // n! but without the p factors
    [[nodiscard]] u64 factorial(u64 n) const
    {
        u64 ret = 1;

        while (n > 0)
        {
            ret = (u64)((__int128)ret * factorialWithoutP(n) % PE);
            n /= P;
        }

        return ret;
    }

    // O(log(n) * e)
    [[nodiscard]] u64 binomial(u64 n, u64 r) const
    {
        if (r < 0 || r > n)
            return 0;
        u64 const binom_padic = factorialValuation(n, P) - factorialValuation(r, P) - factorialValuation(n - r, P);
        if ((i64)binom_padic >= E)
            return (u64)0;

        return (u64)((u128)pPower[binom_padic] * factorial(n) % PE *
                     modInverse((u64)((u128)factorial(r) * factorial(n - r) % PE), PE) % PE * pstirling[binom_padic] %
                     PE);
    }

    u64 operator()(u64 n, u64 r) const { return binomial(n, r); }
};
} // namespace euler
