#include "euler.hpp"
#include "euler/BinomialModPrimePower.hpp"

#include "common.hpp"

using namespace std;
using Int = int64_t;

inline void testWithRandomInputs(auto f, int maxBinaryDigits = 32)
{
    mt19937_64 rng;
    for (int d = 2; d <= maxBinaryDigits; ++d)
    {
        uniform_int_distribution<int64_t> dist(-((int64_t)1 << d), (int64_t)1 << d);
        // cout << "Testing " << d << " binary digits" << endl;
        for (int i = 0; i < 1000; ++i)
            f(rng, dist);
    }
}

inline void testModInverse()
{
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int64_t const a = dist(rng);
            int64_t const m = abs(dist(rng)) + 1;
            if (gcd(a, m) != 1)
                return;
            auto const result = modInverse(a, m);
            auto const product = mod(modmul(result, a, m), m);
            if ((m == 1 && product != 0) || (m != 1 && product != 1))
                fail("modInverse: "s + to_string(a) + " mod " + to_string(m));
        },
        62);
    pass("modInverse");
}

inline void testPowm()
{
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int64_t const a = dist(rng);
            int64_t const b = abs(dist(rng));
            int64_t const c = abs(dist(rng)) + 2;
            auto result = powm(a, b, c);
            auto expected = boost::multiprecision::powm(a, b, c);
            if (result != expected)
                fail("powm");
        },
        31);
    pass("pow");
}

inline void testPowmSafe()
{
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int64_t const a = dist(rng);
            int64_t const b = abs(dist(rng));
            int64_t const c = abs(dist(rng)) + 2;
            auto result = powmSafe(a, b, c);
            auto expected = boost::multiprecision::powm(a, b, c);
            if (result != expected)
                fail("powmSafe");
            return true;
        },
        62);
    pass("powmSafe");
}

// inline void testModPow()
// {
// 	test([](auto &&rng, auto &&dist)
// 		 {
// 			 int64_t a = dist(rng), b = abs(dist(rng)), m = abs(dist(rng));
// 			 if (m == 0)
// 			 	return true;
// 			 auto result = powm(a, b, m);
// 			 auto a2 = (cpp_int)a;
// 			 auto expected = boost::multiprecision::pow(a2, b) % m;
// 			 if (result != expected)
// 			 {
// 				 cout << "Failed: expected " << a << "^" << b << " mod " << m << " to give "
// 				 	  << expected << " but got " << result << endl;
// 				 return false;
// 			 }
// 			 return true; },
// 		 15);
// }

inline void testCrt()
{
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            int64_t const m = abs(dist(rng));
            int64_t const n = abs(dist(rng));
            if (gcd(m, n) != 1 || m == 0 || n == 0)
                return;
            // auto result = crt(array{a, b}, array{m, n});
            auto result = crt(a, b, m, n);
            if (mod(result, m) != mod(a, m) || mod(result, n) != mod(b, n))
                fail("CRT");
        },
        31);
    pass("CRT");
}

inline void testMobius()
{
    int const limit = 1e7;
    SPF const spfs(limit);
    auto mu = mobiusSieve(limit, spfs);
    vector<pair<int, int8_t>> const cases = {{1, 1},
                                             {2, -1},
                                             {3, -1},
                                             {4, 0},
                                             {5, -1},
                                             {6, 1},
                                             {30, -1},
                                             {210, 1},
                                             {100, 0},
                                             {9999998, 1},
                                             {9999999, 0},
                                             {10'000'000, 0},
                                             {1'000'000'000, 0},
                                             {999'999'995, 1}};
    for (auto [n, m] : cases)
        if (n < (int)mu.size() && mu[n] != m)
            fail("mobiusSieve");
    pass("Mobius sieve");
}

inline void testBinom()
{
    std::vector<std::vector<Int>> C(100, std::vector<Int>(100, 0));
    C[0][0] = 1;

    Int mod = (1e7 + 19);
    mod *= mod;

    BinomialModPrimePower const binom(1e7 + 19, 2);

    for (int i = 1; i < 100; ++i)
    {
        C[i][0] = 1;
        for (int j = 1; j <= i; ++j)
        {
            C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % mod;
            assertEqual(C[i][j], binom(i, j));
        }
    }
    pass("Binomial mod prime power");
}

inline void testUtils()
{
    assert(sum(1, 100) == 5050);
    assert(sum(1, (int64_t)100) == 5050);
    assert(sum((int64_t)1, 100) == 5050);
    assert(product(1, 6) == 720);
    assert(product((int64_t)1, 6) == 720);
    assert(product((int64_t)1, (int128_t)6) == 720);
    pass("Sum and product");
}

inline void testEnumCombinations()
{
    Int total = 0;
    it::combinations(range(1, 100), 3)([&](auto &&) { ++total; });
    assert(total == 161700);
    pass("Enum combinations");
}

inline void testBisections()
{
    auto f = [](Int x) { return x * (x + 1) / 2; };
    auto g = [](Int x) { return x / 1000; };

    assert(bisectionLowerBound(f, 15, 0, 1'000'000) == 5);
    assert(bisectionUpperBound(f, 15, 0, 1'000'000) == 6);
    assert(bisectionLowerBound(f, 16, 0, 1'000'000) == 6);
    assert(bisectionUpperBound(f, 16, 0, 1'000'000) == 6);
    assert(bisectionLowerBound(g, 15, 0, 1'000'000) == 15000);
    assert(bisectionUpperBound(g, 15, 0, 1'000'000) == 16000);
    pass("Bisections");
}

inline void testIsPrime()
{
    for (int n = 0; n < 100'000; ++n)
        assert(isPrime(n) == boost::multiprecision::miller_rabin_test(n, 8));
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            auto n = abs(dist(rng));
            if (isPrime(n) != boost::multiprecision::miller_rabin_test(n, 8))
                fail("isPrime");
        },
        62);
    pass("isPrime");
}

inline void testMultiplicativeOrder()
{
    constexpr int N = 1'000'000;
    SPF spfs{N};
    auto totients = totientSieve(N);
    auto result = sum(execution::par, 3, N, [&](int64_t i) {
        removeFactors(i, 2);
        removeFactors(i, 5);
        return i == 1 ? (int64_t)0 : multiplicativeOrder((int64_t)10, i, totients[i], spfs);
    });
    assertEqual(result, 55535191115);
    pass("multiplicativeOrder");
}

inline void testIsSquare()
{
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int64_t const a = abs(dist(rng)) + 2;
            bool const ok = !isSquare(a * a + 1) && !isSquare(a * a - 1) && isSquare(a * a);
            if (!ok)
                fail("isSquare");
        },
        31);
    pass("isSquare");
}

inline void testIsqrt()
{
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            Int const a = abs(dist(rng)) + 1;
            bool const ok = isqrt(a * a + 1) == a && isqrt(a * a - 1) == a - 1 && isqrt(a * a) == a;
            if (!ok)
                fail("isqrt");
        },
        31);
    pass("isqrt");
}

inline void testIt()
{
    assert(it::range(1, 100).sum() == 5050);
    assert(it::range(1, 100).map([](auto &&x) { return x + 1; }).sum() == 5150);
    pass("it::range");
}

inline void testModmul()
{
    static constexpr int m = 1'000'000'007;
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int const a = dist(rng);
            int const b = dist(rng);
            bool const ok = modmul(a, b, m) == (Int)a * b % m;
            if (!ok)
                fail("modmul is wrong for a = "s + to_string(a));
        },
        31);
    pass("modmul");
}

inline void testZMod()
{
    using Z = ZMod<1'000'000'007>;
    testWithRandomInputs(
        [](auto &&rng, auto &&dist) {
            int const a = dist(rng);
            int const b = dist(rng);
            Z const x = Z(a) * b;
            bool const ok = x.value() == mod((Int)a * b, Z::modulus);
            if (!ok)
                fail("ZMod is wrong for a = "s + to_string(a));
        },
        31);
    pass("ZMod");
}

inline void testModInverseUnsigned()
{
    uint64_t const a = 235098327;
    uint64_t const b = 1'000'000'007;
    uint64_t const res = modInverse(a, b);
    assertEqual(res * a % b, 1);
    pass("testModInverseUnsigned");
}

inline void testModUnsignedModulus()
{
    uint32_t const modulus = 61;
    for (int i = -126; i <= 126; ++i)
    {
        auto const res = mod(i, modulus);
        assertEqual(res, (i + 5 * modulus) % modulus);
    }
    pass("testModUnsignedModulus");
}

int main()
{
    auto t1 = now();
    testModUnsignedModulus();
    testModmul();
    testZMod();
    testIt();
    testIsSquare();
    testModInverse();
    testPowm();
    testPowmSafe();
    testCrt();
    testMobius();
    testUtils();
    testEnumCombinations();
    testBisections();
    testIsPrime();
    testMultiplicativeOrder();
    testIsqrt();
    testBinom();
    testModInverseUnsigned();
    println(now() - t1);
}
