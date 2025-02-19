#include "euler.hpp"
#include "euler/BinomialModPrimePower.hpp"

#include "common.hpp"

using namespace std;
using Int = int64_t;

constexpr auto passStr = "\033[0;32;1m[PASS]\033[0m "sv;
constexpr auto failStr = "\033[0;31;1m[FAIL]\033[0m "sv;

inline bool testWithRandomInputs(auto passPred, int maxBinaryDigits = 32)
{
    mt19937_64 rng;
    for (int d = 2; d <= maxBinaryDigits; ++d)
    {
        bool passed = true;
        uniform_int_distribution<int64_t> dist(-((int64_t)1 << d), (int64_t)1 << d);
        // cout << "Testing " << d << " binary digits" << endl;
        for (int i = 0; i < 1000; ++i)
        {
            passed = passPred(rng, dist);
            if (!passed)
                break;
        }
        if (!passed)
        {
            cout << failStr << "First fail at " << d << " binary digits" << '\n';
            return false;
        }
    }
    return true;
}

inline void testModInv()
{
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int64_t a = dist(rng);
                int64_t m = abs(dist(rng)) + 1;
                if (m == 0 || gcd(a, m) != 1)
                    return true;
                auto result = modInverse(a, m);
                auto product = mod(modmul(result, a, m), m);
                if (product != mod((int64_t)1, m))
                {
                    cout << failStr << "modInverse(" << a << ", " << m << ") gave " << result << " but " << result
                         << " * " << a << " = " << product << '\n';
                    return false;
                }
                return true;
            },
            62))
        cout << passStr << "modInverse\n";
}

inline void testPowm()
{
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int64_t a = dist(rng);
                int64_t b = abs(dist(rng));
                int64_t c = abs(dist(rng)) + 2;
                auto result = powm(a, b, c);
                auto expected = boost::multiprecision::powm(a, b, c);
                if (result != expected)
                {
                    cout << failStr << "powm" << tuple{a, b, c} << " gave " << result << " but expected " << expected
                         << "\n";
                    return false;
                }
                return true;
            },
            31))
        cout << passStr << "powm\n";
}

inline void testPowmSafe()
{
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int64_t a = dist(rng);
                int64_t b = abs(dist(rng));
                int64_t c = abs(dist(rng)) + 2;
                auto result = powmSafe(a, b, c);
                auto expected = boost::multiprecision::powm(a, b, c);
                if (result != expected)
                {
                    cout << failStr << "powmSafe" << tuple{a, b, c} << " gave " << result << " but expected "
                         << expected << "\n";
                    return false;
                }
                return true;
            },
            62))
        cout << passStr << "powmSafe\n";
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
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int64_t a = dist(rng);
                int64_t b = dist(rng);
                int64_t m = abs(dist(rng));
                int64_t n = abs(dist(rng));
                if (gcd(m, n) != 1 || m == 0 || n == 0)
                    return true;
                // auto result = crt(array{a, b}, array{m, n});
                auto result = crt(a, b, m, n);
                if (mod(result, m) != mod(a, m) || mod(result, n) != mod(b, n))
                {
                    cout << failStr << "crt(" << a << ", " << b << ", " << m << ", " << n << ") gave " << result
                         << '\n';
                    return false;
                }
                return true;
            },
            31))
        cout << passStr << "CRT\n";
    else
        cout << failStr << "CRT\n";
}

inline void testMobius()
{
    int limit = 1e7;
    SPF spfs(limit);
    auto mu = mobiusSieve(limit, spfs);
    vector<pair<int, int8_t>> cases = {{1, 1},
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
    bool passed = true;
    for (auto [n, m] : cases)
        if (n < (int)mu.size() && mu[n] != m)
        {
            cout << failStr << "at n = " << n << '\n';
            passed = false;
        }
    if (passed)
        cout << passStr << "mobius sieve" << '\n';
}

inline void testBinom()
{
    std::vector<std::vector<Int>> C(100, std::vector<Int>(100, 0));
    C[0][0] = 1;

    Int mod = (1e7 + 19);
    mod *= mod;

    BinomialModPrimePower binom(1e7 + 19, 2);

    for (int i = 1; i < 100; ++i)
    {
        C[i][0] = 1;
        for (int j = 1; j <= i; ++j)
        {
            C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % mod;
            assert(C[i][j] == binom(i, j));
        }
    }
    cout << passStr << "Binomial mod prime power" << '\n';
}

inline void testUtils()
{
    assert(sum(1, 100) == 5050);
    assert(sum(1, (int64_t)100) == 5050);
    assert(sum((int64_t)1, 100) == 5050);
    assert(product(1, 6) == 720);
    assert(product((int64_t)1, 6) == 720);
    assert(product((int64_t)1, (int128_t)6) == 720);
    cout << passStr << "sum and product" << '\n';
}

inline void testEnumCombinations()
{
    Int total = 0;
    it::combinations(range(1, 100), 3)([&](auto &&) { ++total; });
    assert(total == 161700);
    cout << passStr << "Enum combinations\n";
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
    cout << passStr << "Bisections\n";
}

inline void testIsPrime()
{
    for (int n = 0; n < 100'000; ++n)
        assert(isPrime(n) == boost::multiprecision::miller_rabin_test(n, 8));
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                auto n = abs(dist(rng));
                return isPrime(n) == boost::multiprecision::miller_rabin_test(n, 8);
            },
            62))
        cout << passStr << "isPrime\n";
}

inline void testFactor()
{
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                auto n = abs(dist(rng)) + 1;
                return evalPF(factor(n)) == n;
            },
            48))
        cout << passStr << "testFactor\n";
}

inline void testMultiplicativeOrder()
{
    constexpr int N = 1'000'000;
    auto spfs = spfSieve(N);
    auto totients = totientSieve(N);
    auto result = sum(execution::par, 3, N, [&](int64_t i) {
        valuationDivide(i, 2);
        valuationDivide(i, 5);
        return i == 1 ? (int64_t)0 : multiplicativeOrder((int64_t)10, i, totients[i], spfs);
    });
    if (result == 55535191115)
        cout << passStr << "multiplicativeOrder\n";
    else
        cout << failStr << "multiplicativeOrder: " << result << ", expected " << 55535191115 << '\n';
}

inline void testIsSquare()
{
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int64_t a = abs(dist(rng)) + 2;
                bool ok = !isSquare(a * a + 1) && !isSquare(a * a - 1) && isSquare(a * a);
                if (!ok)
                {
                    cout << failStr << "isSquare is wrong for a = " << a << "\n";
                    return false;
                }
                return true;
            },
            31))
        cout << passStr << "isSquare\n";
    else
        cout << failStr << "isSquare\n";
}

inline void testIsqrt()
{
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                Int a = abs(dist(rng)) + 1;
                bool ok = isqrt(a * a + 1) == a && isqrt(a * a - 1) == a - 1 && isqrt(a * a) == a;
                if (!ok)
                {
                    cout << failStr << "isqrt is wrong for a = " << a << "\n";
                    return false;
                }
                return true;
            },
            31))
        cout << passStr << "isqrt\n";
    else
        cout << failStr << "isqrt\n";
}

inline void testIt()
{
    assert(it::range(1, 100).sum() == 5050);
    assert(it::range(1, 100).map([](auto &&x) { return x + 1; }).sum() == 5150);
    cout << passStr << "itertools\n";
}

inline void testModmul()
{
    static constexpr int m = 1'000'000'007;
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int a = dist(rng);
                int b = dist(rng);
                bool ok = modmul(a, b, m) == (Int)a * b % m;
                if (!ok)
                {
                    cout << failStr << "modmul is wrong for a = " << a << "\n";
                    return false;
                }
                return true;
            },
            31))
        cout << passStr << "modmul\n";
    else
        cout << failStr << "modmul\n";
}

inline void testZMod()
{
    using Z = ZMod<1'000'000'007>;
    if (testWithRandomInputs(
            [](auto &&rng, auto &&dist) {
                int a = dist(rng);
                int b = dist(rng);
                Z x = Z(a) * b;
                bool ok = x.value() == mod((Int)a * b, Z::modulus);
                if (!ok)
                {
                    cout << failStr << "ZMod is wrong for a = " << a << "\n";
                    return false;
                }
                return true;
            },
            31))
        cout << passStr << "ZMod\n";
    else
        cout << failStr << "ZMod\n";
}

inline void testFloorsArray()
{
    constexpr Int N = 2e9;
    int L = (int)(0.25 * pow(N, 2.0 / 3));
    assert(floors_array<>::sumTotient(N)[N] == 1215854204348393714);
    for (int sieveSize : {100, 999, 1001, (int)isqrt(N), L - 1, L, L + 1})
    {
        auto sieve = totientSieve(sieveSize);
        assert(floors_array<>::sumTotient(N, sieve)[N] == 1215854204348393714);
    }
    cout << passStr << "floors_array\n";
}

int main()
{
    auto t1 = now();
    testFloorsArray();
    testModmul();
    testZMod();
    testIt();
    testIsSquare();
    testModInv();
    testPowm();
    testPowmSafe();
    testCrt();
    testMobius();
    testUtils();
    testEnumCombinations();
    testBisections();
    testIsPrime();
    testFactor();
    testMultiplicativeOrder();
    testIsqrt();
    testBinom();
    println(now() - t1);
}
