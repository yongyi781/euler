#include "euler/math.hpp"

#include "common.hpp"

using namespace std;

int main()
{
    assertEqual(sum(totientSieve(INT64_C(1))), 1);
    assertEqual(sum(totientSieve(INT64_C(10))), 32);
    assertEqual(sum(totientSieve(INT64_C(100))), 3044);
    assertEqual(sum(totientSieve(INT64_C(1000))), 304192);
    assertEqual(sum(totientSieve(INT64_C(10000))), 30397486);
    assertEqual(sum(totientSieve(INT64_C(100000))), 3039650754);
    assertEqual(sum(totientSieve(INT64_C(1000000))), 303963552392);
    assertEqual(sum(totientSieve(INT64_C(10000000))), 30396356427242);
    assertEqual(sum(totientSieve(INT64_C(100000000))), 3039635516365908);
    pass("totient_sieve");
}
