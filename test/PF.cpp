#include "euler/PF.hpp"

#include "common.hpp"
#include "euler/it/factor.hpp"

using namespace std;

void testPF1()
{
    PF const pf1{{{2, 3}, {3, 7}, {5, 9}}};
    PF const pf2{{{3, 7}, {5, 4}, {7, 11}}};

    PF const expectedMul{{{2, 3}, {3, 14}, {5, 13}, {7, 11}}};
    PF const expectedDiv{{{2, 3}, {5, 5}, {7, -11}}};
    PF const expectedLCM{{{2, 3}, {3, 7}, {5, 9}, {7, 11}}};
    PF const expectedGCD{{{3, 7}, {5, 4}}};

    assert(pf1 != pf2);

    assertEqual(pf1 * pf2, expectedMul);
    assertEqual(pf1 / pf2, expectedDiv);
    assertEqual(pf1 | pf2, expectedLCM);
    assertEqual(pf1 & pf2, expectedGCD);
    assertEqual(pf1.value(), 34171875000);
    assertEqual(pf2.value(), 2702758491838125);
}

void testPF2()
{
    for (int i = 1; i <= 1000; ++i)
    {
        auto const pfi = factor(i);
        assertEqual(pfi.value(), i);
        assertEqual(pfi.totient(), totient(i));
        for (int j = 1; j <= 1000; ++j)
        {
            auto const pfj = factor(j);
            assertEqual(pfi.divides(pfj), j % i == 0);
        }
    }
}

int main()
{
    testPF1();
    testPF2();
    pass("PF");
}
