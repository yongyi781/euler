#include "euler/SegmentTree.hpp"

#include "test.hpp"

using namespace std;

struct SumNode
{
    int value;

    constexpr SumNode(int v = 0) : value(v) {}

    constexpr SumNode operator+(const SumNode &other) const { return {value + other.value}; }

    bool operator<=>(const SumNode &other) const = default;
    friend std::ostream &operator<<(std::ostream &o, const SumNode &n) { return o << n.value; }
};

int main()
{
    // Test case 1: Basic initialization and query
    SegmentTree<SumNode> const st1(array{1, 2, 3, 4, 5});
    assertEqual(st1.size(), 5);
    assertEqual(st1.query(0, 5).value, 15);
    assertEqual(st1.query(0, 3).value, 6);
    assertEqual(st1.query(2, 5).value, 12);

    // Test case 2: Update and query
    SegmentTree<SumNode> st2(array{1, 2, 3, 4, 5});
    st2.update(0, SumNode(10));
    assertEqual(st2.query(0, 5).value, 24);
    st2.update(4, SumNode(10));
    assertEqual(st2.query(0, 5).value, 29);

    // Test case 3: Empty tree
    SegmentTree<SumNode> const st3(0);
    assertEqual(st3.size(), 0);

    // Test case 4: Single element tree
    SegmentTree<SumNode> st4(array{100});
    assertEqual(st4.size(), 1);
    assertEqual(st4.query(0, 1).value, 100);
    st4.update(0, SumNode(200));
    assertEqual(st4.query(0, 1).value, 200);

    // Test case 5: All method
    SegmentTree<SumNode> const st5(array{1, 2, 3, 4, 5});
    assertEqual(st5.all().value, 15);

    // Test case 6: operator[]
    SegmentTree<SumNode> const st6(array{1, 2, 3});
    assertEqual(st6[0].value, 1);
    assertEqual(st6[1].value, 2);
    assertEqual(st6[2].value, 3);

    // Test case 7: Query with i == j
    SegmentTree<SumNode> const st7(array{1, 2, 3});
    assertEqual(st7.query(0, 0).value, 0); // Querying an empty range should return identity element

    // Test case 8: Query with i > j (should be handled by query logic, typically returns identity)
    // This case is not explicitly handled by the current query, but for a monoid,
    // an empty range should result in the identity element. SumNode's default constructor
    // provides this (value = 0).
    assertEqual(st7.query(2, 1).value, 0);

    // Test case 9: Large number of elements
    std::vector<SumNode> large_data;
    large_data.reserve(1000);
    for (int i = 0; i < 1000; ++i)
        large_data.emplace_back(i + 1);
    SegmentTree<SumNode> const st8(large_data);
    assertEqual(st8.size(), 1000);
    assertEqual(st8.query(0, 1000).value, 500500); // Sum of 1 to 1000

    // Test case 10: Update multiple times
    SegmentTree<SumNode> st9(array{1, 1, 1, 1, 1});
    st9.update(0, SumNode(2));
    st9.update(1, SumNode(3));
    st9.update(2, SumNode(4));
    assertEqual(st9.query(0, 5).value, 11); // 2+3+4+1+1 = 11

    // Test case 11: Query with non-zero start index
    SegmentTree<SumNode> const st10(array{10, 20, 30, 40, 50});
    assertEqual(st10.query(1, 3).value, 50); // 20 + 30

    // Test case 12: Query with range spanning across internal nodes
    SegmentTree<SumNode> const st11(array{1, 2, 3, 4, 5, 6, 7, 8});
    assertEqual(st11.query(1, 7).value, 27); // 2+3+4+5+6+7
}
