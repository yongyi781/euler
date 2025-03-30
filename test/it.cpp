#include "euler/it.hpp"

#include "common.hpp"
#include "euler/ZMod.hpp"

using namespace std;

void testRange()
{
    assertEqual(it::range(1, 100).size(), 100);
    assertEqual(it::range(1, 100, 3).size(), 34);
    assertEqual(it::range(1, 99, 3).size(), 33);
    assertEqual(it::range(1, 100).map([&](auto &&x) { return x / 2; }).unique().size(), 51);
    pass("it/range");
}

void testCombinatorics()
{
    assertEqual(it::combinations(vector{1, 2, 3}, 2).size(), 3);
    assertEqual(it::combinations(it::range(1, 30).to(), 15).size(), 155117520);
    pass("it/combinatorics");
}

void testDigits()
{
    assertEqual(it::digits(123).to(), vector{3, 2, 1});
    assertEqual(it::digits(123, 10).to(), vector{3, 2, 1});
    assertEqual(it::digits(123, 2).to(), vector{1, 1, 0, 1, 1, 1, 1});
    pass("it/digits");
}

void testDFinite()
{
    using Z = ZMod<(int64_t)1'234'567'891>;
    assertEqual(
        827177246,
        it::dfinite(
            {{160, 464, 480, 208, 32}, {-88, -364, -432, -200, -32}, {-108, -181, -110, -27, -2}, {36, 81, 62, 19, 2}},
            vector<Z>{0, 2, 6}, 1)[10000]);
    pass("it/dfinite");
}

int main()
{
    testRange();
    testCombinatorics();
    testDigits();
    testDFinite();
    pass("it");
}
