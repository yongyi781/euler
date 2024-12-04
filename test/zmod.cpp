#include "euler/ZMod.hpp"

#include "common.hpp"

using namespace std;

int main()
{
    assertEqual(ZMod<11>(3).pow(2), 9);
    pass("ZMod");
}
