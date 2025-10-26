#include "euler/literals.hpp"
#include "euler/math.hpp"

#include "test.hpp"

using namespace std;

void testFloorDiv()
{
    i64 const N = 100;
    for (i64 a = 1; a <= N; ++a)
    {
        for (i64 b = 1; b <= N; ++b)
        {
            if (b == 0)
                continue;
            assertEqual(floorDiv(a, b), a / b, tuple{"floor", a, b});
            assertEqual(ceilDiv(a, b), (a + b - 1) / b, tuple{"ceil", a, b});

            assertEqual(floorDiv(-a, -b), a / b, tuple{"floor", -a, -b});
            assertEqual(ceilDiv(-a, -b), (a + b - 1) / b, tuple{"ceil", -a, -b});

            assertEqual(floorDiv(-a, b), (-a - mod(-a, b)) / b, tuple{"floor", -a, b});
            assertEqual(ceilDiv(-a, b), (-a + mod(a, b)) / b, tuple{"ceil", -a, b});

            assertEqual(floorDiv(a, -b), (-a - mod(-a, b)) / b, tuple{"floor", a, -b});
            assertEqual(ceilDiv(a, -b), (-a + mod(a, b)) / b, tuple{"ceil", a, -b});
        }
    }

    // Now random testing
    mt19937_64 rng;
    uniform_int_distribution dist_a(1_i64, 1_i64 << 55);
    uniform_int_distribution dist_b(1_i64, 1_i64 << 8);
    i64 const trials = 10000;
    for (i64 i = 0; i < trials; ++i)
    {
        i64 const a = dist_a(rng);
        i64 b = dist_b(rng);
        while (b == 0)
            b = dist_b(rng);
        assertEqual(floorDiv(a, b), a / b, tuple{"floor", a, b});
        assertEqual(ceilDiv(a, b), (a + b - 1) / b, tuple{"ceil", a, b});

        assertEqual(floorDiv(-a, -b), a / b, tuple{"floor", -a, -b});
        assertEqual(ceilDiv(-a, -b), (a + b - 1) / b, tuple{"ceil", -a, -b});

        assertEqual(floorDiv(-a, b), (-a - mod(-a, b)) / b, tuple{"floor", -a, b});
        assertEqual(ceilDiv(-a, b), (-a + mod(a, b)) / b, tuple{"ceil", -a, b});

        assertEqual(floorDiv(a, -b), (-a - mod(-a, b)) / b, tuple{"floor", a, -b});
        assertEqual(ceilDiv(a, -b), (-a + mod(a, b)) / b, tuple{"ceil", a, -b});
    }
    pass("floorDiv and ceilDiv");
}

void testFastDiv()
{
    i64 const N = 100;
    for (i64 a = 1; a <= N; ++a)
    {
        for (i64 b = 1; b <= N; ++b)
        {
            if (b == 0)
                continue;
            assertEqual(fastDiv(a, b), a / b, tuple{"fastDiv", a, b});
            assertEqual(fastDiv(-a, -b), a / b, tuple{"fastDiv", -a, -b});
            assertEqual(fastDiv(-a, b), (-a - mod(-a, b)) / b, tuple{"fastDiv", -a, b});
            assertEqual(fastDiv(a, -b), (-a - mod(-a, b)) / b, tuple{"fastDiv", a, -b});
        }
    }

    // Now random testing
    mt19937_64 rng;
    uniform_int_distribution dist_a(1_i64, 1_i64 << 55);
    uniform_int_distribution dist_b(1_i64, 1_i64 << 8);
    i64 const trials = 10000;
    for (i64 i = 0; i < trials; ++i)
    {
        i64 const a = dist_a(rng);
        i64 b = dist_b(rng);
        while (b == 0)
            b = dist_b(rng);
        assertEqual(fastDiv(a, b), a / b, tuple{"fastDiv", a, b});
        assertEqual(fastDiv(-a, -b), a / b, tuple{"fastDiv", -a, -b});
        assertEqual(fastDiv(-a, b), (-a - mod(-a, b)) / b, tuple{"fastDiv", -a, b});
        assertEqual(fastDiv(a, -b), (-a - mod(-a, b)) / b, tuple{"fastDiv", a, -b});
    }
    pass("fastDiv");
}

int main()
{
    testFloorDiv();
    testFastDiv();
    pass("floor_div");
}
