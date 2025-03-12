#include "euler/PF.hpp"

#include "common.hpp"

using namespace std;

int main()
{
    PF const pf1{{{2, 3}, {3, 7}, {5, 9}}};
    PF const pf2{{{3, 7}, {5, 4}, {7, 11}}};

    PF const expectedMul{{{2, 3}, {3, 14}, {5, 13}, {7, 11}}};
    PF const expectedDiv{{{2, 3}, {5, 5}, {7, -11}}};
    PF const expectedLCM{{{2, 3}, {3, 7}, {5, 9}, {7, 11}}};
    PF const expectedGCD{{{3, 7}, {5, 4}}};

    assert(pf1 != pf2);

    // assertEqual(pf1 * pf2, expectedMul);
    assertEqual(pf1 / pf2, expectedDiv);
    assertEqual(pf1 | pf2, expectedLCM);
    assertEqual(pf1 & pf2, expectedGCD);
    pass("dirichlet");
}
