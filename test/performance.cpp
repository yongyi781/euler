#include <algorithm>
#include <iostream>

#include "euler.hpp"
#include <boost/rational.hpp>

using namespace std;
using Int = int64_t;

/** Safe but slower version. */
template <typename T> constexpr auto powmSafe2(const T &base, const integral2 auto &exponent, integral2 auto modulus)
{
    assert(modulus > 0);
    if constexpr (integral2<T>)
    {
        return exponent < 0 ? modInverse(pow(base % modulus, -exponent, T(1), mod_multiplies_safe(modulus)), modulus)
                            : pow(base % modulus, exponent, T(1), mod_multiplies_safe(modulus));
    }
    else
    {
        assert(exponent >= 0);
        return pow(base % modulus, exponent, T(1), mod_multiplies_safe(modulus));
    }
}

template <typename F1, typename F2> void compare(F1 f1, F2 f2, int trials)
{
    mt19937_64 rng;
    uniform_int_distribution<int64_t> dist(-(int64_t)1 << 22, (int64_t)1 << 22);
    int128_t total = 0;
    auto t1 = now();
    for (int i = 0; i < trials; ++i)
        total += (int128_t)f1(rng, dist);
    auto t2 = now();
    cout << "Trials = " << trials << " | Method 1 = " << (t2 - t1) / trials << " | ans = " << total << '\n';

    rng = mt19937_64();
    total = 0;
    t1 = now();
    for (int i = 0; i < trials; ++i)
        total += (int128_t)f2(rng, dist);
    t2 = now();
    cout << "Trials = " << trials << " | Method 2 = " << (t2 - t1) / trials << " | ans = " << total << '\n';
}

// 586 ns and 283 ns with clang++, 390 ns and 288 ns with g++.
inline void comparePowm(int trials = 1000000)
{
    cout << "Comparing boost powm and my powm\n";
    compare(
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = abs(dist(rng)) + 1;
            auto modulus = abs(dist(rng)) + 1;
            return boost::multiprecision::powm(abs(x), y, modulus);
        },
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = abs(dist(rng)) + 1;
            auto modulus = abs(dist(rng)) + 1;
            return powm(abs(x), y, modulus);
        },
        trials);
}

// 586 ns and 283 ns with clang++, 390 ns and 288 ns with g++.
inline void comparePowmSafe(int trials = 1000000)
{
    cout << "Comparing boost powm and my powmSafe\n";
    compare(
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = abs(dist(rng)) + 1;
            auto modulus = abs(dist(rng)) + 1;
            return boost::multiprecision::powm(abs(x), y, modulus);
        },
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = abs(dist(rng)) + 1;
            auto modulus = abs(dist(rng)) + 1;
            return powmSafe(abs(x), y, modulus);
        },
        trials);
}

// 22 ns for both (g++), 31 ns for both (clang++).
inline void comparePowInt128(int trials = 10000000)
{
    cout << "Comparing boost pow and my pow on int128s\n";
    compare(
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = mod(dist(rng), 20);
            return boost::multiprecision::pow((boost::multiprecision::int128_t)x, y);
        },
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = mod(dist(rng), 20);
            return pow((int128_t)x, y);
        },
        trials);
}

// 237 ns and 254 ns with clang++, 286 ns and 349 ns with g++.
// New crt only takes 162 ns now. Woo! 228 ns with g++. (6/22/2024).
inline void compareCrt(int trials = 1000000)
{
    cout << "Copmaring 4-argument CRT and 2-argument CRT\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            int64_t const m = abs(dist(rng)) + 1;
            int64_t const n = abs(dist(rng)) + 1;
            if (gcd(m, n) != 1 || m == 0 || n == 0)
                return (int64_t)0;
            auto result = crt(a, b, m, n);
            return result;
        },
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            int64_t const m = abs(dist(rng)) + 1;
            int64_t const n = abs(dist(rng)) + 1;
            if (gcd(m, n) != 1 || m == 0 || n == 0)
                return (int64_t)0;
            auto result = crt(array{a, b}, array{m, n});
            return result;
        },
        trials);
}

// 6 ns and 25 ns with clang++, 5 ns and 29 ns with g++.
inline void compareSqrt(int trials = 10000000)
{
    cout << "Copmaring boost sqrt and isqrt\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = abs(dist(rng));
            return boost::multiprecision::sqrt(a);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = abs(dist(rng));
            return isqrt(a);
        },
        trials);
}

// 6 ns and 25 ns with clang++, 5 ns and 29 ns with g++.
inline void compareSqrt2(int trials = 10000000)
{
    cout << "Copmaring sqrt + cast and isqrt\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = abs(dist(rng));
            return (int64_t)sqrt(a);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = abs(dist(rng));
            return isqrt(a);
        },
        trials);
}

// 27 ns and 57 ns with clang++, 25 ns and 50 ns with g++.
inline void compareCbrt(int trials = 10000000)
{
    cout << "Copmaring boost cbrt and std::cbrt\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            return (int64_t)boost::math::cbrt(a);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            return (int64_t)cbrt(a);
        },
        trials);
}

// 47 ns and 64 ns with clang++, 63 ns and 90 ns with g++.
inline void compareGcd(int trials = 10000000)
{
    cout << "Copmaring std::gcd and boost::integer::gcd\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            return std::gcd(a, b);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            return boost::integer::gcd(a, b);
        },
        trials);
}

// 51 ms and 68 ns with clang++, 72 ns and 97 ns with g++.
inline void compareLcm(int trials = 10000000)
{
    cout << "Copmaring std::lcm and boost::integer::lcm\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            return std::lcm(a, b);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = dist(rng);
            return boost::integer::lcm(a, b);
        },
        trials);
}

template <int N> void compareConstantPow(int trials = 10000000)
{
    cout << "Copmaring std::pow and boost::integer::pow\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t a = dist(rng);
            return (int64_t)pow((double)a, N);
        },
        [](auto &rng, auto &dist) {
            int64_t a = dist(rng);
            return (int64_t)boost::math::pow<N>((double)a);
        },
        trials);
}

auto floorDiv2(integral2 auto a, integral2 auto b)
{
    using T = decltype(a);
    T const d = a / b;
    T const r = a % b;
    return r ? (d - ((a < 0) ^ (b < 0))) : d;
}

void compareFloorDiv(int trials = 10000000)
{
    cout << "Copmaring floorDiv and floorDiv2\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = abs(dist(rng)) + 1;
            return floorDiv(a, b);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = dist(rng);
            int64_t const b = abs(dist(rng)) + 1;
            return floorDiv2(a, b);
        },
        trials);
}

// powm faster with primes, modInverse faster with prime powers.
inline void compareModInv(int trials = 10000000)
{
    cout << "Copmaring modInverse and powm\n";
    compare(
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            Int const modulus = static_cast<Int>(10627) * 10627;
            if (x % modulus == 0)
                return (int64_t)0;
            return modInverse(x, modulus);
        },
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            Int const modulus = static_cast<Int>(10627) * 10627;
            if (x % modulus == 0)
                return (int64_t)0;
            return powm(x, modulus * (modulus - 1) - 1, modulus);
        },
        trials);
}

inline void compareInt128(int trials = 10000000)
{
    cout << "Copmaring __int128 and boost::multiprecision::int128_t\n";
    compare(
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = abs(dist(rng));
            auto z = (__int128)x * 9223372036854775808ULL + y;
            return z % 324235235;
        },
        [](auto &rng, auto &dist) {
            auto x = dist(rng);
            auto y = abs(dist(rng));
            auto z = (int128_t)x * 9223372036854775808ULL + y;
            return z % 324235235;
        },
        trials);
}

// powm faster with primes, modInverse faster with prime powers.
inline void compareFastDiv(int trials = 10000000)
{
    cout << "Copmaring integer division and fastDiv\n";
    compare(
        [](auto &rng, auto &dist) {
            int64_t const a = abs(dist(rng));
            int64_t const b = abs(dist(rng)) + 1;
            return floorDiv(a, b);
        },
        [](auto &rng, auto &dist) {
            int64_t const a = abs(dist(rng));
            int64_t const b = abs(dist(rng)) + 1;
            return fastDiv(a, b);
        },
        trials);
}

// Building an array of 30k factorials took:
// mpz_int: 282 ms
// cpp_int: 1.45 s
// sagemath: 101 ms??? Fastest? How?
// Python.: 600 ms
// VC++ cpp_int: 1.67 s.
template <typename T> auto factorialTest1()
{
    constexpr int N = 30000;
    vector<T> fact(N + 1);
    fact[0] = 1;

    for (int i = 1; i <= N; ++i)
        fact[i] = fact[i - 1] * i;

    return fact[15];
}

// cpp_int: 757 ms (g++), 734 ms (clang++)
// int128_t: 175 ms (clang++)
template <typename T> auto factorialTest2()
{
    vector<T> output(30000000);
    output[0] = 1;

    transform_inclusive_scan(execution::par, counting_iterator(1), counting_iterator((int)output.size()),
                             output.begin() + 1, mod_multiplies(1000000007), [](int i) { return T(i); });
    // auto output = unfold(30000000, T(1), mod_multiplies(1000000007));
    return output[29999999];
}

/**
 * __int128: 1.46 s (g++ and clang++)
 * int128_t: 1.90 s (g++ and clang++)
 * int256_t: 386.808 ms (g++), 595 ms (clang++)
 * int512_t: 372.812 ms (g++), 601 ms (clang++)
 * cpp_int: 499.97 ms (g++), 1.85 s (clang++)
 * VC++ int128_t: 1.0 s
 * VC++ cpp_int: 2.0 s
 */
template <typename T> auto performanceTest2()
{
    constexpr int N = 100'000'000;
    constexpr int MOD = 1'000'000'007;
    T num = 1;

    for (int64_t i = 1; i < N; ++i)
    {
        num *= (2 * i + 1);
        num %= MOD;
    }

    return num;
}

inline auto performanceTest3()
{
    constexpr int N = 10'000'000;
    mt19937_64 rng;
    vector<int256_t> v(N);

    generate_n(execution::par, v.begin(), N, [&] { return rng(); });
    return v[100];
}

inline auto performanceTest4()
{
    constexpr int64_t N = 100'000'000;
    auto v = randomVector(N, 0, 10, mt19937_64{random_device{}()});
    return v[100];
}

// Simple for loop multiplying the first 10^8 cubes
// int128_t: 141 ms
// int256_t: 1.2 s
// int512_t: 1.6 s
// int1024_t: 1.8 s
// cpp_int: 2.06 s with g++, 2.65 s with clang++
// mpz_int: 12.7 s
// Sage math (Integer): 2.74 s for 10^7
// Sage math (int): 774 ms for 10^7
// Python: 900 ms for 10^7
// Numpy (int): 1.4 s for 10^7
inline auto performanceTest5()
{
    boost::multiprecision::cpp_int total = 0;
    for (boost::multiprecision::cpp_int i = 0; i < 100'000'000; ++i)
        total += i * i * i;
    return total;
}

template <typename T> std::vector<T> spf2(T limit)
{
    std::vector<T> sieve(limit + 1);
    for (T i = 2; i <= limit; ++i)
        sieve[i] = i % 2 == 0 ? 2 : i % 3 == 0 ? 3 : i % 5 == 0 ? 5 : i;
    T sqrtLimit = isqrt(limit);
    if (limit >= 2000000)
    {
        for (T i = 7; i <= sqrtLimit; ++i)
            if (sieve[i] == i)
                for_each(std::execution::par, counting_iterator(i), counting_iterator(limit / i) + 1, [&sieve, i](T j) {
                    if (sieve[i * j] == i * j)
                        sieve[i * j] = i;
                });
    }
    else
    {
        for (T i = 7; i <= sqrtLimit; ++i)
            if (sieve[i] == i)
                for (T j = i * i; j <= limit; j += i)
                    if (sieve[j] == j)
                        sieve[j] = i;
    }
    return sieve;
}

// 468 ms with g++ and clang++.
inline auto spfTest()
{
    auto v = spfSieve((int)1e8);
    return v[17051];
}

inline auto spfTest2()
{
    auto v = spf2((int)1e8);
    return v[17051];
}

inline auto powmTest1()
{
    int4096_t const base = 2;
    int4096_t const exponent = ((int4096_t)1) << 4095;
    int4096_t const modulus = (int64_t)1e18;
    auto result = boost::multiprecision::powm(base, exponent, modulus);
    return result;
}

inline auto powmTest2()
{
    int128_t const base = 2;
    int4096_t const exponent = ((int4096_t)1) << 4095;
    // cout << exponent << endl;
    int64_t const modulus = 1e18;
    auto result = ::powm(base, exponent, modulus);
    return result;
}

inline auto multiplyPFTest()
{
    Factorization<int> const pf1{};
    Factorization<int> const pf2{{2, 4}, {7, 4}, {11, 6}, {13, 8}};
    auto a = mulPF(pf1, pf2);
}

vector<int> spfs;

inline void factorTest()
{
    cout << "Doing the factor test\n";
    spfs = spfSieve(10000000);
    for (int i = 1; i < 10000000; ++i)
        factor(i, spfs);
}

inline void digitsTest()
{
    int256_t const n("1234567");
    auto s = to_string(n);
    vector<int> digits(s.size());
    std::ranges::transform(s, digits.begin(), [](char c) { return (int)c - '0'; });
    std::ranges::reverse(digits);
}

int main()
{
    comparePowm();
    comparePowmSafe();
    comparePowInt128();
    compareCrt();
    compareSqrt();
    compareSqrt2();
    compareFloorDiv();
    compareFastDiv();
    printTiming(spfTest);
    printTiming(spfTest2);
    printTiming(powmTest1);
    printTiming(powmTest2);
    printTiming(multiplyPFTest);
    printTiming(factorTest);
    printTiming(digitsTest);
}
