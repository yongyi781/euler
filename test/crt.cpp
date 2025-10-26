#include "euler/math.hpp"

#include "test.hpp"

using namespace std;

int main()
{
    assertEqual(crt(9, 24, 332, 123), 31881);
    assertEqual(crt(0, 5, 1, 7), 5);
    assertEqual(crt(vector{1}, vector{5}), 1);
    assertEqual(crt(vector{1, 2, 3}, vector{7, 8, 9}), 498);
}
