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
/*  if(argc < 6) { return 1; }
  double bias = std::stod(std::string(argv[1]));
  std::cout << "bias: " << bias << '\n';
  double delta = std::stod(std::string(argv[2]));
  std::cout << "delta: " << delta << '\n';
  uint32_t toConsider = std::stoul(std::string(argv[3]));
  std::cout << "toConsider: " << toConsider << '\n';
  double shepardExponent = std::stod(std::string(argv[4]));
  std::cout << "shepExp: " << shepardExponent << '\n';
  double avRelSize = std::stod(std::string(argv[5]));
  std::cout << "avRelSize: " << avRelSize << '\n';*/

/*using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 1u, Data, 4, 3>;
  typename ShepIntpol::DataTransfer data;
  double number = 0.0;
  double const max = 100.0;
  for(double number = 0.0; number < max; number += delta) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    // auto sum = number * number;          // relerr<0.007 if bias: 4321 delta: 2.3 toConsider: 4 shepExp: 3 TL: 3
    auto sum = std::sin(number / 7);        // biaserr<0.05 if bias: 2 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.006 unbiaserr<0.06 if bias: 12 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.003 unbiaserr<0.033 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.002 unbiaserr<0.04 if bias: 42 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.002 unbiaserr<0.033 if bias: 22 delta: 1.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.002 unbiaserr<0.033 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 2
                                            // biaserr<0.0007 unbiaser<0.007 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0007 unbiaser<0.007 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 4
                                            // optimal bias = |vertspan| * 22
// biaserr is what seen on biased relative error graph
    item.mPayload = {static_cast<double>(sum)};
    data.push_back(item);
  }
  ShepIntpol shep(data, toConsider, avRelSize, shepardExponent, bias);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  double d2 = 0.25;
  std::cout << "diff=[";
  for(double i = 0u; i < max; i+=d2) {
    typename ShepIntpol::Location loc{i};
    // auto v = i*i;
    auto v = std::sin(i / 7);
    valSum += v;
    auto d = (v == 0.0 ? 1.0 : shep.interpolate(loc)[0] / v);
    diffSum += d * d;
    std::cout << d << (i < max - d2 ? ", " : "];\n");
  }
  auto diffAvg = std::sqrt(diffSum) / max / max;
  auto valAvg = valSum / max / max;
//  std::cout << "err: " << diffAvg/valAvg << '\n';
*/
  /*using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 2u, Data, 6, 3>;
                                            // biaserr<0.0007 unbiaser<0.015 if bias: 17000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0008 unbiaser<0.015 if bias: 15000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0006 unbiaser<0.015 if bias: 19000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0006 unbiaser<0.015 if bias: 19000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 2
                                            // optimal bias = |vertspan| * 22
  typename ShepIntpol::DataTransfer data;
  double const max = 64.0;
  for(double n1 = 0.0; n1 < max; n1 += delta) {
    for(double n2 = 0.0; n2 < max; n2 += delta) {
      typename ShepIntpol::Data item;
      item.mLocation[0] = n1;
      item.mLocation[1] = n2 * 10.0;
      auto sum = (n1 * n1 + n2 * n2) / 10.0;
      item.mPayload = {sum};
      data.push_back(item);
    }
  }
  ShepIntpol shep(data, toConsider, avRelSize, shepardExponent, bias);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  std::cout << "diff=[";
  uint32_t limit = 64u;
  for(uint32_t i = 0u; i < limit; ++i) {
    for(uint32_t j = 0u; j < limit; ++j) {
      typename ShepIntpol::Location loc{i, j * 10.0};
      auto v = (i*i+j*j)/10.0;
      valSum += v;
      auto d = (v == 0.0 ? 1.0 : shep.interpolate(loc)[0] / v);
      diffSum += d * d;
      std::cout << d << (j < limit - 1u ? ", " : (i < limit - 1u ? ";\n" : "];\n"));
    }
  }
  auto diffAvg = std::sqrt(diffSum) / limit / limit;
  auto valAvg = valSum / limit / limit;
  std::cout << "err: " << diffAvg/valAvg << '\n';*/

/*  double slice = std::stod(std::string(argv[6]));
  std::cout << "slice: " << avRelSize << '\n';
  using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 3u, Data, 6, 3>;

  typename ShepIntpol::DataTransfer data;
  double const max = 64.0;
  for(double n1 = 0.0; n1 < max; n1 += delta) {
    for(double n2 = 0.0; n2 < max; n2 += delta) {
      for(double n3 = 0.0; n3 < max; n3 += delta) {
        typename ShepIntpol::Data item;
        item.mLocation[0] = n1;
        item.mLocation[1] = n2 *  5.0;
        item.mLocation[2] = n3 * 20.0;
        auto sum = (n1 * n1 + n2 * n2 + n3 * n3) / 10.0;
        item.mPayload = {sum};
        data.push_back(item);
      }
    }
  }
  ShepIntpol shep(data, toConsider, avRelSize, shepardExponent, bias);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  std::cout << "diff=[";
  uint32_t limit = 64u;
  for(uint32_t i = 0u; i < limit; ++i) {
    for(uint32_t j = 0u; j < limit; ++j) {
      typename ShepIntpol::Location loc{slice, i * 5.0, j * 20.0};
      auto v = (slice*slice+i*i+j*j)/10.0;
      valSum += v;
      auto d = (v == 0.0 ? 1.0 : shep.interpolate(loc)[0] / v);
      diffSum += d * d;
      std::cout << d << (j < limit - 1u ? ", " : (i < limit - 1u ? ";\n" : "];\n"));
    }
  }
  auto diffAvg = std::sqrt(diffSum) / limit / limit;
  auto valAvg = valSum / limit / limit;
  std::cout << "err: " << diffAvg/valAvg << '\n';
*/
  ShepardRayBending t28(28.0);

  std::vector<double> angles;
  for(auto i = 89.6; i <= 89.95; i += 0.025) {
    angles.push_back(i);
  }

  // TODO check angles
  test("height1dist50height", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 50).first; });
  test("height1dist50dir", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 50).second; });
  test("height1dist100height", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 100).first; });
  test("height1dist100dir", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 100).second; });
  test("height1dist200height", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 200).first; });
  test("height1dist200dir", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 200).second; });
  test("height1dist300height", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 300).first; });
  test("height1dist300dir", angles, [&t28](auto angle){ return t28.getHeightDirectionBending(1.0, (angle - 90.0) * cgPi / 180.0, 300).second; });

  test("getBendingHorizDisp", angles, [&t28](auto angle){ return t28.getBendingHorizDisp((angle - 90) * cgPi / 180.0); });
  test("getCriticalDirection", std::vector<double>{{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}}, [&t28](auto height){ return t28.getCriticalDirection(height); });
  test("getAsphaltHitHorizDisp", std::vector<double>{{-0.785, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.017}}, [&t28](auto dir){ return t28.getAsphaltHitHorizDisp(dir); });

  // TODO double getAsphaltHitHorizDisp(double const aDirection) {

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
