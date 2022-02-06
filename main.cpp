#include "Angle2apparentMirrorDepth.h"
#include <functional>
#include <array>
#include <iostream>


/*void test(std::vector<double> const& aSamplesX, uint32_t const aDegreeMinus, uint32_t const aDegreePlus, std::function<double(double)> aF) {
  std::vector<double> y;
  for(auto x : aSamplesX) {
    std::cout << x << ", ";
    y.push_back(aF(x));
  }
  std::cout << '\n';
  for(auto x : aSamplesX) {
    std::cout << aF(x) << ", ";
  }
  std::cout << '\n';
  PolynomApprox(aSamplesX, y, aDegreeMinus, aDegreePlus);
  std::cout << '\n';
}*/

int main(int argc, char **argv) {
  /*test({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}, 1u, 2u, [](double const aX){ return aX; });
  test({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}, 1u, 2u, [](double const aX){ return aX * aX; });
  test({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}, 2u, 1u, [](double const aX){ return 1.0 / aX; });
  test({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}, 2u, 0u, [](double const aX){ return 1.0 / aX / aX + 3.0; });
  test({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2}, 2u, 2u, [](double const aX){ return aX * aX+ 1.0 / aX / aX; });
  test({0.0005, 0.01, 0.03, 0.13, 0.55, 0.89}, 2u, 1u, [](double const aX){
    return aX < 0.008 ? 63 : (aX < 0.02 ? 42 : (aX < 0.06 ? 30 : (aX < 0.3 ? 25 : (aX < 0.8 ? 21 : 20))));
  });*/

/*  Angle2apparentMirrorDepth t32(12.0);
  Angle2apparentMirrorDepth t40(20.0);*/
  Angle2apparentMirrorDepth t28(28.0);
/*  Angle2apparentMirrorDepth t56(36.0);
  Angle2apparentMirrorDepth t64(44.0);*/
  std::array<double, 15> h = {{0.0005, 0.01, 0.03, 0.05, 0.1, 0.13, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.89}};
std::cout << "h: ";
  for(auto f : h) {
std::cout << f << ", ";
  }
std::cout << "\nt: ";
  for(auto f : h) {
/*    std::cout << " t32 " << t32.getTempAtHeight(f);
    std::cout << " t40 " << t40.getTempAtHeight(f);*/
    std::cout << t28.getTempAtHeight(f) << ", ";
/*    std::cout << " t56 " << t56.getTempAtHeight(f);
    std::cout << " t64 " << t64.getTempAtHeight(f);*/
  }
std::cout << "\n";
  return 0;
}
