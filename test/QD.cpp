#include "euler/QD.hpp"

#include "common.hpp"

using namespace std;

inline auto test(QD q)
{
    i64 const xmin = 0, xmax = 1'000'000, ymin = 0, ymax = 1'000'000;
    auto const res = q.solve(xmin, xmax, ymin, ymax);

    vector<Vector<i64, 2>> expected;
    for (i64 x = xmin; x <= xmax; ++x)
    {
        i64 const d = (q.e + q.b * x) * (q.e + q.b * x) - 4 * q.c * (q.a * x * x + q.d * x + q.f);
        if (d < 0)
            continue;
        i64 const sqrt_d = isqrt(d);
        if (sqrt_d * sqrt_d != d)
            continue;
        int const sign_c = q.c > 0 ? 1 : -1;
        if (d == 0)
        {
            i64 const yc = -q.b * x - q.e;
            if (sign_c * yc >= ymin && yc % (2 * q.c) == 0 && yc / (2 * q.c) <= ymax)
                expected.emplace_back(x, yc / (2 * q.c));
        }
        else
            for (i64 const yc : {-q.b * x - q.e + sqrt_d, -q.b * x - q.e - sqrt_d})
                if (sign_c * yc >= ymin && yc % (2 * q.c) == 0 && yc / (2 * q.c) <= ymax)
                    expected.emplace_back(x, yc / (2 * q.c));
    }

    cout << q;
    if (res != expected)
    {
        cout << ansi::red << " FAIL" << ansi::reset << '\n';
        cout << setw(12) << "actual = " << res << '\n';
        cout << setw(12) << "expected = " << expected << '\n';
    }
    else
        cout << ansi::green << " PASS" << ansi::reset << " (" << res.size() << ")\n";
}

int main()
{
    test({5, 0, -1, 0, 0, 1});
    test({5, 2, -1, 2, 0, 1});
    test({5, 2, -1, 2, 3, -1});
    test({3, 5, -4, 1, 2, -12});
    test({42, 62, 21, 0, 0, -585});

    test({3, 13, 5, -11, -7, -92});
    for (i64 n = -20; n <= 20; ++n)
        test({5, -2, -1, 0, 0, -n});
    for (i64 n = -20; n <= 20; ++n)
        test({5, -2, -1, 1, 0, -n});
    for (i64 n = -20; n <= 20; ++n)
        test({5, -2, -1, 1, -1, -n});
}
