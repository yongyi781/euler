#pragma once

#include "prime.hpp"

inline namespace euler
{
/// A class to compute binomial coefficients mod prime powers.
class BinomialModPrimePower
{
  public:
    // O(p*min(e, p) + e^2)
    BinomialModPrimePower(int64_t p, int e) : _p(p), _e(e), fact(e * 2 - 1), ifact(e * 2 - 1), pfact(e * 2 - 1)
    {
        pPower.resize(e);

        for (int i = 0; i < e; ++i)
        {
            pPower[i] = _pe;
            _pe *= p;
        }

        int ep = (int)(e < p ? e : p);
        stirling = std::vector<std::vector<int64_t>>(p + 1, std::vector<int64_t>(ep + 1, 0));
        stirling[0][0] = 1;
        for (int i = 1; i <= p; ++i)
            for (int j = 1; j <= ep; ++j)
                stirling[i][j] = (int64_t)(((int128_t)(i - 1) * stirling[i - 1][j] + stirling[i - 1][j - 1]) % _pe);

        prods.resize(e * 2 - 1);
        int64_t prod = 1;
        int64_t invStirling = modInverse(stirling[p][1], _pe);

        for (int i = 0; i <= e * 2 - 3; ++i)
        {
            prods[i] = prod;
            prod = (int64_t)((int128_t)prod * risingFactorial(i, p - 1) % _pe * invStirling % _pe);
        }
        prods[e * 2 - 2] = prod;

        pstirling.resize(e);
        pstirling[0] = 1;
        for (int i = 1; i < e; ++i)
            pstirling[i] = (int64_t)((int128_t)pstirling[i - 1] * stirling[p][1] % _pe);

        int len = e * 2 - 1;
        fact[0] = 1LL;
        pfact[0] = 0;
        for (int i = 1; i < len; ++i)
        {
            int64_t num = i;
            pfact[i] = 0;
            while (num % p == 0)
            {
                ++pfact[i];
                num /= p;
            }

            ifact[i - 1] = num;
            pfact[i] += pfact[i - 1];
            fact[i] = (int64_t)((__int128)fact[i - 1] * num % _pe);
        }

        ifact[len - 1] = modInverse(fact[len - 1], _pe);
        for (int i = len - 2; i >= 0; --i)
            ifact[i] = (int64_t)((__int128)ifact[i + 1] * ifact[i] % _pe);

        lagrangeCoeff.resize(len);
        for (int i = 0; i < len; ++i)
        {
            auto denominator = (int64_t)((__int128)ifact[i] * ifact[len - 1 - i] % _pe);
            if (((len - 1 - i) & 1) != 0)
                denominator = _pe - denominator;

            lagrangeCoeff[i] = (int64_t)((__int128)denominator * prods[i] % _pe);
        }
    }

    [[nodiscard]] constexpr int64_t p() const { return _p; }

    [[nodiscard]] constexpr int64_t e() const { return _e; }

    [[nodiscard]] constexpr int64_t pe() const { return _pe; }

    // n! % p^e, but p, 2p, ... is not multiplied
    [[nodiscard]] int64_t factorialWithoutP(int64_t n) const
    {
        int64_t ndp = n / _p;
        int64_t ret = risingFactorial(ndp, n % _p);

        if (ndp <= _e * 2 - 2)
            return (int64_t)((__int128)ret * prods[ndp] % _pe);

        ret = (int64_t)((__int128)ret * lagrangeInterpolate(ndp) % _pe);
        return ret;
    }

    // n! but without the p factors
    [[nodiscard]] int64_t factorial(int64_t n) const
    {
        int64_t ret = 1;

        while (n > 0)
        {
            ret = (int64_t)((__int128)ret * factorialWithoutP(n) % _pe);
            n /= _p;
        }

        return ret;
    }

    // O(log(n) * e)
    [[nodiscard]] int64_t binomial(int64_t n, int64_t r) const
    {
        if (r < 0 || r > n)
            return 0;
        int64_t binom_padic = factorialValuation(n, _p) - factorialValuation(r, _p) - factorialValuation(n - r, _p);
        if (binom_padic >= _e)
            return 0LL;

        return (int64_t)((int128_t)pPower[binom_padic] * factorial(n) % _pe *
                         modInverse((int64_t)((int128_t)factorial(r) * factorial(n - r) % _pe), _pe) % _pe *
                         pstirling[binom_padic] % _pe);
    }

    int64_t operator()(int64_t n, int64_t r) const { return binomial(n, r); }

  private:
    int64_t _p = 1, _pe = 1;
    int _e = 0;

    // p^k
    std::vector<int64_t> pPower;

    // First kind
    std::vector<std::vector<int64_t>> stirling;

    // (p-1)!^k
    std::vector<int64_t> pstirling;

    // factorials and its inverse, but without the factor p
    std::vector<int64_t> fact, ifact;
    // p-adic valuation of factorial
    std::vector<int> pfact;

    std::vector<int64_t> lagrangeCoeff;

    std::vector<int64_t> prods;

    // Computes (np+1)(np+2) ... (np+m) % p^e
    // with m < p
    [[nodiscard]] int64_t risingFactorial(int64_t n, int64_t m) const
    {
        int64_t ret = 0;
        int64_t pn = 1;
        int ep = (int)(_e < _p ? _e : _p);
        for (int j = 0; j < ep; ++j)
        {
            ret = (int64_t)(((__int128)stirling[m + 1][j + 1] * pn + ret) % _pe);
            pn = (int64_t)((__int128)pn * _p * n % _pe);
        }

        return ret;
    }

    [[nodiscard]] int64_t lagrangeInterpolate(int64_t ndp) const
    {
        int len = _e * 2 - 1;
        int64_t ptot = 0;

        std::vector pfactorsnum(len, 0);
        std::vector prenum(len, 0LL);
        std::vector sufnum(len, 0LL);

        for (int i = 0; i < len; ++i)
        {
            int64_t num = ndp - i;
            pfactorsnum[i] = 0;
            while (num % _p == 0)
            {
                num /= _p;
                ++pfactorsnum[i];
            }
            ptot += pfactorsnum[i];
            prenum[i] = sufnum[i] = num;

            if (i > 0)
                prenum[i] = (int64_t)((__int128)prenum[i - 1] * prenum[i] % _pe);
        }

        for (int i = len - 2; i >= 0; --i)
            sufnum[i] = (int64_t)((__int128)sufnum[i] * sufnum[i + 1] % _pe);

        int64_t sum = 0;
        for (int j = 0; j < len; ++j)
        {
            int j2 = len - 1 - j;
            int64_t pfactor = ptot - pfactorsnum[j] - pfact[j] - pfact[j2];

            if (pfactor >= _e)
                continue;

            auto numerator = (int64_t)((__int128)(j > 0 ? prenum[j - 1] : 1) * (j < len - 1 ? sufnum[j + 1] : 1) % _pe);

            sum = (int64_t)(((__int128)numerator * lagrangeCoeff[j] % _pe * pPower[pfactor] + sum) % _pe);
        }

        return sum;
    }
};
} // namespace euler
