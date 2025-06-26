#include "euler/algorithm.hpp"

#include "test.hpp"

using namespace std;

int main()
{
    vector const v{0, 0, 5, 2, 4, 3, 2, 1, 5, 6};
    vector const expected{0, 0, 5, 7, 11, 14, 16, 17, 22, 28};
    assertEqual(partialSum(execution::par, v), expected);
    assertEqual(partialSum(v), expected);
    pass("partial_sum");
}
