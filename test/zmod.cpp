#include "euler/ZMod.hpp"

#include "test.hpp"

using namespace std;
using namespace euler;

int main() { assertEqual(ZMod<11>(3).pow(2), 9); }
