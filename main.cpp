#include "PolynomialRayBending.h"
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
/*  PolynomialRayBending t28(28.0);
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
  PolynomialRayBending::Gather g;
  g.mAsphalt = false;
  g.mCollection.emplace_back(1.0, 2.0, 2.0/::sqrt(0.1*0.1+2.0*2.0));
  g.mCollection.emplace_back(0.9, 3.0, 3.0/::sqrt(0.1*0.1+3.0*3.0));
  g.mCollection.emplace_back(0.8, 4.0, 4.0/::sqrt(0.1*0.1+4.0*4.0));
  g.mCollection.emplace_back(0.75, 0.0, std::nan("1"));

  auto p = PolynomialRayBending::toRayPath(g);
  std::cout << "d = [";
  for(auto i : p.mCollection) {
    std::cout << i.mHorizDisp << ", ";
  }
  std::cout << "\nh = [";
  for(auto i : p.mCollection) {
    std::cout << i.mHeight << ", ";
  }
  std::cout << "\na = [";
  for(auto i : p.mCollection) {
    std::cout << i.mAngleFromHoriz << ", ";
  }
  std::cout << '\n';
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
