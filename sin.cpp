#include "OdeSolverGsl.h"
#include <iostream>
#include <utility>
#include <cmath>


//clang++ -std=c++17 sin.cpp -o sin -ggdb -lgsl
 

class Sin final {
public:
  static constexpr uint32_t csNvar = 1u; 
  using Real                       = double;
  using Variables                  = std::array<Real, csNvar>; 

  int differentials(double const aX, const double aY[], double aDydt[]) const {
    aDydt[0] = std::cos(aX);
    return GSL_SUCCESS;
  }

  int jacobian(double, const double aY[], double *aDfdy, double aDfdt[]) const { 
    return GSL_SUCCESS;
  }
};

class Sin2 final { // y''=-y  y''+y=0  D=-4 r=+-j  y=cos(x)+jsin(x) 
// dy/dt = z     0 1
// dz/dt = -y
public:
  static constexpr uint32_t csNvar = 2u; 
  using Real                       = double;
  using Variables                  = std::array<Real, csNvar>; 

  int differentials(double const aX, const double aY[], double aDydt[]) const {
    aDydt[0] = aY[1];
    aDydt[1] = -aX;
    return GSL_SUCCESS;
  }

  int jacobian(double, const double aY[], double *aDfdy, double aDfdt[]) const { 
    return GSL_SUCCESS;
  }
};

using Diff = Sin;

int main() {
  const double atol=1.0e-6, rtol=atol, h1=0.01, hmin=0.0, x1=0.0, x2=3.1415926539 * 60, dx = 3.1415926539 / 5.0;
  Diff diffEq;
  OdeSolverGsl solver(StepperType::cRungeKutta23, x1, x2, atol, rtol, h1, diffEq);
  typename Diff::Variables ystart;
  ystart[0]=0.0;
  std::vector<double> ys;
  std::vector<double> xs;
  std::cout << "x=[";
  for (double x = x1; x < x2; x += dx) {
    auto [t, solution] = solver.solve(ystart, [x](double const aX, typename Diff::Variables const&){ return x <= aX; });
    xs.push_back(t);
    ys.push_back(solution[0]);
    std::cout << t << ", ";
  }
  std::cout << "\ny=[";
  for(uint32_t i = 0u; i < ys.size(); ++i) {
    std::cout << ys[i] - std::sin(xs[i]) << ", ";
  }
  std::cout << "\n";
}
