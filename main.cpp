#include "Angle2apparentMirrorDepth.h"
#include <iostream>

#include <functional>

void test(double const aStart, double const aStep, double const aEnd, uint32_t const aDegreeMinus, uint32_t const aDegreePlus, std::function<double(double)> aF) {
  std::vector<double> x;
  std::vector<double> y;
  for(auto i = aStart; i <= aEnd; i += aStep) {
    std::cout << i << ", ";
    x.push_back(i);
    y.push_back(aF(i));
  }
  std::cout << '\n';
  for(auto i : y) {
    std::cout << i << ", ";
  }
  std::cout << '\n';
  PolynomApprox(x, y, aDegreeMinus, aDegreePlus);
  std::cout << '\n';
}

void test(std::vector<double> const &aSampleX, uint32_t const aDegreeMinus, uint32_t const aDegreePlus, std::function<double(double)> aF) {
  std::vector<double> y;
  for(auto x : aSampleX) {
    std::cout << x << ", ";
    y.push_back(aF(x));
  }
  std::cout << '\n';
  for(auto i : y) {
    std::cout << i << ", ";
  }
  std::cout << '\n';
  PolynomApprox(aSampleX, y, aDegreeMinus, aDegreePlus);
  std::cout << '\n';
}

int main(int argc, char **argv) {
/*  test(1.0, 1.0, 7.0, 1u, 2u, [](double const aX){ return aX; });
  test(1.0, 1.0,  7.0, 1u, 2u, [](double const aX){ return aX * aX; });
  test(0.1, 0.1, 0.7, 2u, 1u, [](double const aX){ return 1.0 / aX; });
  test(0.1, 0.02, 0.7, 2u, 1u, [](double const aX){ return 1.0 / aX; });
  test(0.1, 0.1, 0.7, 2u, 0u, [](double const aX){ return 1.0 / aX / aX + 3.0; });
  test(0.1, 0.02, 0.7, 2u, 0u, [](double const aX){ return 1.0 / aX / aX + 3.0; });
  test(0.1, 0.1, 1.9, 2u, 2u, [](double const aX){ return aX * aX + 1.0 / aX / aX; });
  test(0.1, 0.02, 1.9, 2u, 2u, [](double const aX){ return aX * aX + 1.0 / aX / aX; });*/
  test({0.0005, 0.01, 0.03, 0.13, 0.55, 0.89}, 5u, 0u, [](double const aX){
    return aX < 0.008 ? 63 : (aX < 0.02 ? 42 : (aX < 0.06 ? 30 : (aX < 0.3 ? 25 : (aX < 0.8 ? 21 : 20))));
  });

/*  Angle2apparentMirrorDepth t32(20.0, 12.0);
  Angle2apparentMirrorDepth t40(20.0, 20.0);
  Angle2apparentMirrorDepth t48(20.0, 28.0);
  Angle2apparentMirrorDepth t56(20.0, 36.0);*/
  Angle2apparentMirrorDepth t64(20.0, 44.0);
  std::array<double, 12> h = {{0.0005, 0.01, 0.03, 0.07, 0.13, 0.25, 0.38, 0.55, 0.65, 0.77, 0.89, 1.0}};
  for(auto f : h) {
    std::cout << "h " << f;
/*    std::cout << " t32 " << t32.getTempAtHeight(f);
    std::cout << " t40 " << t40.getTempAtHeight(f);
    std::cout << " t48 " << t48.getTempAtHeight(f);
    std::cout << " t56 " << t56.getTempAtHeight(f);*/
    std::cout << " t64 " << t64.getTempAtHeight(f);
    std::cout << '\n';
  }
  return 0;
}
