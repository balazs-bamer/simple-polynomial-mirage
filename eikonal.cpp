#include "nr3.h"
#include "stepper.h"
#include "stepperdopr853.h"
#include "odeint.h"

double refract(double const aH) {
  auto celsius = 24.0 + 0.018 +  6.37 * std::exp(-aH / 10.08); // Toled vettem.
  return 1.0 + 7.86e-4 * 101 / (celsius + 273);
}

double refractDiff(double const aH) {
  return 5.6866e-7 * std::exp(0.0992063 * aH) / std::pow(0.0214465 + std::exp(0.0992063 * aH), 2.0);
}


class Eikonal final {
public:
  static constexpr double cgC = 299792458; // m/2
  
private:
  double mDeltaZ;

public:
  Eikonal(double const aDeltaZ) : mDeltaZ(aDeltaZ) {}

  void operator() (const double/* aS*/, VecDoub_I const &aY, VecDoub_O &aDyds) {
    double n    = refract(aY[2]);
    double v    = cgC / n;
//    double dvdz = -v / n * refractDiff(aY[2]);

    aDyds[0] = v * aY[3];
    aDyds[1] = v * aY[4];
    aDyds[2] = v * aY[5];
    aDyds[3] = 0.0;
    aDyds[4] = 0.0;
    aDyds[5] = refractDiff(aY[2]) / v / n;
    //aDyds[5] = -1.0 / v / v * dvdz;
  }
};

void comp(char * aPre, double const aDirDegree, double const aHeight) {
  const Int nvar = 6;
  const double atol=1.0e-6, rtol=atol, h1=0.01, hmin=0.0, s1=0.0, s2=10000;
  VecDoub ystart(nvar);
  ystart[0] = 0.0;
  ystart[1] = 0.0;
  ystart[2] = aHeight;
  auto u = refract(aHeight) / Eikonal::cgC;
  ystart[3] = u * std::cos(aDirDegree / 180.0 * 3.1415926539);
  ystart[4] = 0.0;
  ystart[5] = u * std::sin(aDirDegree / 180.0 * 3.1415926539);
  Output out(100);
  Eikonal eikonal(0.001);
  Odeint<StepperDopr853<Eikonal>> ode(ystart, s1, s2, atol, rtol, h1, hmin, out, eikonal);
  ode.integrate();
  std::cout << aPre << "x=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.ysave[0][i] << (i < out.count - 1 ? ", " : "];\n");
  }
  std::cout << aPre << "z=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.ysave[2][i] << (i < out.count - 1 ? ", " : "];\n");
  }
  std::cout << "\n";
}

int main() {
  comp("i_005_1", -0.05, 1.0);
  comp("i_002_1", -0.02, 1.0);
  comp("i_001_1", -0.01, 1.0);
}
