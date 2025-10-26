#include "euler/find_recurrence.hpp"
#include "euler/io.hpp"
#include <boost/rational.hpp>

#include "test.hpp"

using namespace std;

int main()
{
    vector<int64_t> const v{3, 10, 23, 44, 77, 128, 205, 320, 491, 744, 1117, 1666, 2473, 3658, 5397};
    vector<int64_t> const expected{1, -2, 2, -3, 3, -1};
    assertEqual(findRecurrence(v, (int64_t)1'000'000'007), expected);

    vector<boost::rational<int64_t>> v2{3, 10, 23, 44, 77, 128, 205, 320, 491, 744, 1117, 1666, 2473, 3658, 5397};
    vector<boost::rational<int64_t>> const expected2{1, -2, 2, -3, 3, -1};
    assertEqual(findRecurrence(v2), expected2);
}
