#include "PolynomialRayBending.h"
#include "simpleRaytracer.h"
#include <functional>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>


constexpr double cgEpsilon = 0.0001f;

bool eq(double const aF1, double const aF2, double const aEpsilon) {
  return ::abs(aF1 - aF2) < aEpsilon;
}

bool eq(double const aF1, double const aF2) {
  return eq(aF1, aF2, cgEpsilon);
}

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

  std::vector<double> xs, ys, zs;
  double const deltay = 0.25;
  double const deltax = 50;
  double const r2    = 25.0;
  std::cout << "tx = [";
  for(double x = -800; x < 1190; x += deltax) {
    std::cout << x << ", ";
  }
  std::cout << "1200]';\n";
  std::cout << "ty = [";
  for(double x = -8.0; x < 1.9; x += deltay) {
    std::cout << x << ", ";
  }
  std::cout << "2]';\n";
  std::cout << "ref = [";
  for(double x = -800; x < 1210; x += deltax) {
    for(double y = -8.0; y < 2.1; y += deltay) {
      double ex = (x - 200) / 200;
      double ey = y + 3;
      double e = ::exp(ex / 10.0 + ey / 10.0) * 2.0;
      ex *= ex; ey *= ey;
      double z;
      if(ex + ey < r2) {
        z = ::sqrt(r2 - ex - ey) + e;
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
      }
      else z = 0.0;
      std::cout << z;
      if(eq(y, 2)) {
        if(eq(x, 1200))
          std::cout << "];\n";
        else
          std::cout << ";\n";
      }
      else
        std::cout << ',';
    }
  }

  PolynomApprox poly(zs, {{xs.data(), 5u}, {ys.data(), 5u}});
  std::cout << "err = " << poly.getRrmsError() << "\n";
  std::cout << "ev = [";
  for(double x = -800; x < 1210; x += deltax) {
    for(double y = -8.0; y < 2.1; y += deltay) {
      double ex = (x - 200) / 200;
      double ey = y + 3;
      double e = ::exp(ex / 10.0 + ey / 10.0) * 2.0;
      ex *= ex; ey *= ey;
      double z;
      if(ex + ey < r2) {
        z = ::sqrt(r2 - ex - ey) + e - poly.eval(std::initializer_list<double>{x ,y});
      }
      else z = 0.0;
      std::cout << z;
      if(eq(y, 2)) {
        if(eq(x, 1200))
          std::cout << "];\n";
        else
          std::cout << ";\n";
      }
      else
        std::cout << ',';
    }
  }

/*  PolynomialRayBending t28(28.0);

  std::vector<double> angles({89.6, 89.65, 89.7, 89.75, 89.8, 89.85, 89.9, 89.95});

  test("height1dist50height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 50).first; });
  test("height1dist50dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 50).second; });
  test("height1dist100height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 100).first; });
  test("height1dist100dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 100).second; });
  test("height1dist200height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 200).first; });
  test("height1dist200dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 200).second; });
  test("height1dist300height", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 300).first; });
  test("height1dist300dir", angles, [&t28](auto angle){ return t28.getHeightDirection(1.0, (90.0 - angle) * cgPi / 180.0, 300).second; });*/
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
