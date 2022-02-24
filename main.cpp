#include "Angle2apparentMirrorDepth.h"
#include "simpleRaytracer.h"
#include <functional>
#include <array>
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
  Angle2apparentMirrorDepth t28(28.0);
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
  test("refractionAtHeight", heights, [&t28](auto aHeight){ return t28.getRefractionAtHeight(aHeight); });
  test("refractionAtTemp", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0}, [&t28](auto aTemp){ return t28.getRefractionAtTempRise(aTemp); });
  test("horizDisp", inclinations, [&t28](auto aIncl){ return t28.approximateHorizontalDisplacement(aIncl); });
  test("reflectionDepth", inclinations, [&t28](auto aIncl){ return t28.approximateReflectionDepth(aIncl); });
  test("iterations", inclinations, [&t28](auto aIncl){ return t28.approximateIterations(aIncl); });

  int32_t const cQuantiles = 10;
  auto layerThicknesses = t28.getQuantiles(cQuantiles);
  std::cout << "Layer thickness quantiles (" << cQuantiles << "): ";
  for(auto t : layerThicknesses) {
    std::cout << t << ", ";
  }
  std::cout << "\naverage thickness: " << t28.getAverage() << '\n';

  // TODO parse args
  // Object(char const * const aName, double const aCenterX, double const aCenterY, double const aWidth, double const aHeight)
  Object object(argv[1], 1000.0, 1.8, 2.0, 1.5); // Top is 2.3

  //Medium(double const aTempDiff, Object const &aObject) : mHotPlate(aTempDiff), mObject(aObject) {}
  Medium medium(28.0, object);

  /*Image::Image(double const aCenterX, double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium const &aMedium)*/
//  Image image(0.0, 1.1, 0.0, 0.03*1000.0/(1.0+3.0+1.8), 0.03/40, 80, 40, 1, medium);
  Image image(0.0, 1.05, 0.0, 0.03*1000.0/8, 0.03/8000, 3000, 8000, 4, medium);
  image.process(argv[2]);
  return 0;
}
