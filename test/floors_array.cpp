#include "euler/math.hpp"

#include "common.hpp"

using namespace std;

int main()
{
    for (int64_t n = 1; n <= 10000; ++n)
    {
        auto val = floors_array<>::sumTotient(n)[n];
        auto val2 = sum(1, n, [&](int64_t k) { return totient(k); });
        assertEqual(val, val2);
    }
    pass("floors_array");
}
