#include "euler/dirichlet.hpp"

#include "test.hpp"

using namespace std;

int main()
{
    int64_t const N = 200000;
    auto phi = totientSieve(N);
    vector<int64_t> phiSum(N + 1);
    inclusive_scan(phi.begin(), phi.end(), phiSum.begin());
    for (int64_t n = 1; n <= N; ++n)
    {
        auto expected = phiSum[n];
        auto actual = dirichlet::totient(n).back();
        assertEqual(actual, expected);
    }

    int64_t const N2 = pow((int64_t)10, 9);
    assertEqual(dirichlet::totient(N2)[N2], 303963551173008414);
    pass("dirichlet");
}
