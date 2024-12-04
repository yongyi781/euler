#include "euler/it.hpp"

#include "common.hpp"

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

int main()
{
    testRange();
    testCombinatorics();
    testDigits();
    pass("it");
}
