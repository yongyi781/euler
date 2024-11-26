#include "euler/math.hpp"

#include "common.hpp"

using namespace std;

int main()
{
    assertEqual(crt(9, 24, 332, 123), 31881);
    assertEqual(crt(0, 5, 1, 7), 5);
    assertEqual(crt(vector{1LL}, vector{5LL}), 1LL);
    assertEqual(crt(vector{1LL, 2LL, 3LL}, vector{7LL, 8LL, 9LL}), 498LL);
    pass("CRT");
}
