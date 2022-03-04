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

  std::vector<double> xs, ys, zs, as;
  double const deltaX = 50;
  double const deltaY = 0.25;
  double const deltaZ = 0.125;
  double const minX = -800.0;
  double const midX =  200.0;
  double const maxX = 1200.0;
  double const minY = -8.0;
  double const midY = -3.0;
  double const maxY = 2.0;
  double const minZ = -5.0;
  double const maxZ = 5.0;
  double const r    = 5.0;
  double const r2    = r * r;
  double const zSlice = -2.0;
  std::cout << "tx = [";
  for(double x = minX; x < maxX - cgEpsilon; x += deltaX) {
    std::cout << x << ", ";
  }
  std::cout << maxX << "]';\n";
  std::cout << "ty = [";
  for(double x = minY; x < maxY - cgEpsilon; x += deltaY) {
    std::cout << x << ", ";
  }
  std::cout << maxY << "]';\n";
  std::cout << "tz = [";
  for(double x = minZ; x < maxZ - cgEpsilon; x += deltaZ) {
    std::cout << x << ", ";
  }
  std::cout << maxZ << "]';\n";
  std::cout << "ref = [";
  for(double x = minX; x < maxX + cgEpsilon; x += deltaX) {
    for(double y = minY; y < maxY + cgEpsilon; y += deltaY) {
      for(double z = minZ; z < maxZ + cgEpsilon; z += deltaZ) {
        double ex = (x - midX) / 200;
        double ey = y - midY;
        double e = ::exp(ex / 10.0 - ey / 10.0 + z / 10.0);
        ex *= ex; ey *= ey;
        double ez = z * z;
        double a;
        if(ex + ey + ez <= r2) {
          a = ::sqrt(r2 - ex - ey - ez) + e;
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
          as.push_back(a);
        }
        else a = 0.0;
        if(eq(z, zSlice)) {
          std::cout << a;
          if(eq(y, maxY)) {
            if(eq(x, maxX))
              std::cout << "];\n";
            else
              std::cout << ";\n";
          }
          else
            std::cout << ',';
        }
      }
    }
  }

  PolynomApprox poly(as, {{xs.data(), 2u}, {ys.data(), 2u}, {zs.data(), 2u}});
  std::cout << "err = " << poly.getRrmsError() << "\n";
  std::cout << "ev = [";
  for(double x = minX; x < maxX + cgEpsilon; x += deltaX) {
    for(double y = minY; y < maxY + cgEpsilon; y += deltaY) {
      for(double z = minZ; z < maxZ + cgEpsilon; z += deltaZ) {
        double ex = (x - midX) / 200;
        double ey = y - midY;
        double e = ::exp(ex / 10.0 - ey / 10.0 + z / 10.0);
        ex *= ex; ey *= ey;
        double ez = z * z;
        double a;
        if(ex + ey + ez <= r2) {
          a = ::sqrt(r2 - ex - ey - ez) + e - poly.eval(std::initializer_list<double>{x, y, z});
        }
        else a = 0.0;
        if(eq(z, zSlice)) {
          std::cout << a;
          if(eq(y, maxY)) {
            if(eq(x, maxX))
               std::cout << "];\n";
            else
              std::cout << ";\n";
          }
          else
            std::cout << ',';
        }
      }
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
