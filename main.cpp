#include "ShepardRayBending.h"
#include "simpleRaytracer.h"
#include <functional>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>


void test(char const * const aName, std::vector<double> const& aSamplesX, std::function<double(double)> aF) {
  std::cout << aName << "_x = [";
  for(auto x : aSamplesX) {
    std::cout << std::setprecision(11) << std::scientific << x << (x == aSamplesX.back() ? "];\n" : ", ");
  }
  std::cout << aName << "_y = [";
  for(auto x : aSamplesX) {
    std::cout << std::setprecision(11) << std::scientific << aF(x) << (x == aSamplesX.back() ? "];\n" : ", ");
  }
}

int main(int argc, char **argv) {
  using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 2u, Data, 3>;
  std::vector<ShepIntpol::Data> data;
  uint32_t number = 0u;
  for(uint32_t i = 0u; i < 60u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    auto sum = number * number;
    number = (number + 81u) % 127u;
    item.mLocation[1] = number;
    sum += number * number;
    number = (number + 81u) % 127u;
    item.mPayload = {static_cast<double>(sum)};
    data.push_back(item);
  }
  ShepIntpol shep(data, 3u, argc);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  for(uint32_t i = 0u; i < 127u; ++i) {
    for(uint32_t j = 0u; j < 127u; ++j) {
      typename ShepIntpol::Location loc{i, j};
      auto v = i*i+j*j;
      valSum += v;
      auto d = v - shep.interpolate(loc)[0];
      diffSum += d * d;
//    std::cout << "dftc: " << shep.getDistanceFromTargetCenter(loc) << " diffInt: " << d << '\n';
    }
  }
  auto diffAvg = ::sqrt(diffSum) / 127.0 / 127.0;
  auto valAvg = valSum / 127.0 / 127.0;
  std::cout << "err: " << diffAvg/valAvg << '\n';

  // TODO Octave output of interpolated function.

  /*ShepardRayBending t28(28.0);

  std::vector<double> angles({89.6, 89.65, 89.7, 89.75, 89.8, 89.85, 89.9, 89.95});

  test("height1dist50height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 50).first; });
  test("height1dist50dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 50).second; });
  test("height1dist100height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 100).first; });
  test("height1dist100dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 100).second; });
  test("height1dist200height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 200).first; });
  test("height1dist200dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 200).second; });
  test("height1dist300height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 300).first; });
  test("height1dist300dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 300).second; });
*/  
/*
  std::vector<double> heights({0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.0035, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.13, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8, 0.89});
  std::vector<double> inclinations;
  int const count = 97;
  double actual = t28.getCriticalInclination();
  double incr = (cgPi / 2.0 - actual) / (count - 1);
  for(auto i = 0; i < count; ++i) {
    inclinations.push_back(actual);
    actual += incr;
  }
  test("tempAtHeight", heights, [&t28](auto aHeight){ return t28.getTempRiseAtHeight(aHeight); });
*/
  /*
  // TODO parse args
  // Object(char const * const aName, double const aCenterX, double const aCenterY, double const aWidth, double const aHeight)
  Object object(argv[1], 1000.0, 1.8, 2.0, 1.5); // Top is 2.3

  //Medium(double const aTempDiff, Object const &aObject) : mHotPlate(aTempDiff), mObject(aObject) {}
  Medium medium(28.0, object);

  //Image::Image(double const aCenterX, double const aCenterY,
  //    double const aTilt, double const aPinholeDist,
  //    double const aPixelSize,
  //    uint32_t const aResZ, uint32_t const aResY,
  //    uint32_t const aSubSample, Medium const &aMedium)
//  Image image(0.0, 1.1, 0.0, 0.03*1000.0/(1.0+3.0+1.8), 0.03/40, 80, 40, 1, medium);
  Image image(0.0, 1.05, 0.0, 0.03*1000.0/8, 0.03/8000, 3000, 8000, 4, medium);
  image.process(argv[2]);*/
  return 0;
}
