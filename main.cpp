#include "Angle2apparentMirrorDepth.h"
#include <functional>
#include <array>
#include <iostream>


void test(char const * const aName, std::vector<double> const& aSamplesX, std::function<double(double)> aF) {
  std::cout << aName << "x = [";
  for(auto x : aSamplesX) {
    std::cout << std::scientific << x << ", ";
  }
  std::cout << '\n' << aName << "y = [";
  for(auto x : aSamplesX) {
    std::cout << std::scientific << aF(x) << ", ";
  }
  std::cout << '\n';
}

int main(int argc, char **argv) {
  Angle2apparentMirrorDepth t28(28.0);
  std::vector<double> heights({0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.0035, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05, 0.1, 0.13, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.89});
  test("tempAtHeight", heights, [&t28](auto aHeight){ return t28.getTempAtHeight(aHeight); });
  test("refractionAtHeight", heights, [&t28](auto aHeight){ return t28.getRefractionAtHeight(aHeight); });
  test("refractionAtTemp", {297.5, 300, 305, 310, 315, 320, 325}, [&t28](auto aTemp){ return t28.getRefractionAtTemp(aTemp); });
  return 0;
}
