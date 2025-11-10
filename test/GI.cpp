#include "euler/GI.hpp"
#include "euler/prime.hpp"
#include "test.hpp"

using namespace std;
using namespace euler;

void testGIDivision()
{
    // Test 1: Basic division with no remainder
    GI<int> const g1_t1{5, 0};
    GI<int> const g2_t1{1, 0};
    auto [q1_t1, r1_t1] = g1_t1.div(g2_t1);
    assertEqual(GI<int>{5, 0}, q1_t1, "Test 1 Quotient");
    assertEqual(GI<int>{0, 0}, r1_t1, "Test 1 Remainder");

    // Test 2: Division with a remainder
    GI<int> const g1_t2{5, 5};
    GI<int> const g2_t2{2, 0};
    auto [q1_t2, r1_t2] = g1_t2.div(g2_t2);
    assertEqual(GI<int>{2, 2}, q1_t2, "Test 2 Quotient");
    assertEqual(GI<int>{1, 1}, r1_t2, "Test 2 Remainder");

    // Test 3: Division by a purely imaginary Gaussian integer
    GI<int> const g1_t3{10, 0};
    GI<int> const g2_t3{0, 2};
    auto [q1_t3, r1_t3] = g1_t3.div(g2_t3);
    assertEqual(GI<int>{0, -5}, q1_t3, "Test 3 Quotient");
    assertEqual(GI<int>{0, 0}, r1_t3, "Test 3 Remainder");

    // Test 4: Division where the quotient is zero
    GI<int> const g1_t4{1, 1};
    GI<int> const g2_t4{10, 0};
    auto [q1_t4, r1_t4] = g1_t4.div(g2_t4);
    assertEqual(GI<int>{0, 0}, q1_t4, "Test 4 Quotient");
    assertEqual(GI<int>{1, 1}, r1_t4, "Test 4 Remainder");

    // Test 5: Division by one
    GI<int> const g1_t5{7, -3};
    GI<int> const g2_t5{1, 0};
    auto [q1_t5, r1_t5] = g1_t5.div(g2_t5);
    assertEqual(GI<int>{7, -3}, q1_t5, "Test 5 Quotient");
    assertEqual(GI<int>{0, 0}, r1_t5, "Test 5 Remainder");

    // Test 6: Division by itself
    GI<int> const g1_t6{4, 2};
    GI<int> const g2_t6{4, 2};
    auto [q1_t6, r1_t6] = g1_t6.div(g2_t6);
    assertEqual(GI<int>{1, 0}, q1_t6, "Test 6 Quotient");
    assertEqual(GI<int>{0, 0}, r1_t6, "Test 6 Remainder");

    // Test 7: Division by zero (expecting an exception)
    GI<int> const g1_t7{1, 1};
    GI<int> const g2_t7{0, 0};
    bool caught_exception = false;
    try
    {
        auto _ = g1_t7.div(g2_t7);
    }
    catch (const std::overflow_error &e)
    {
        caught_exception = true;
    }
    assertEqual(true, caught_exception, "Test 7 Exception Handling");

    // Original test case
    GI<int> const g1_orig{-6, 0};
    GI<int> const g2_orig{1, 2};
    auto [q_orig, r_orig] = g1_orig.div(g2_orig);
    assertEqual(GI<int>{-1, 2}, q_orig, "Original Test Quotient");
    assertEqual(GI<int>{-1, 0}, r_orig, "Original Test Remainder");
}

void testPrimeNormGI()
{
    for (auto const p : primeRange(100000))
    {
        if (p % 4 == 3)
            continue;
        auto const g = primeNormGI(p);
        assertEqual(g.norm(), p, p);
    }
}

int main()
{
    testGIDivision();
    testPrimeNormGI();
}
