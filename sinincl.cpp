#include "nr3.h"
#include "stepper.h"
#include "stepperdopr853.h"
#include "odeint.h"

double n(double const aH) {
  auto celsius = 24.0 + 0.018 +  6.37 * std::exp(-aH / 10.08); // Toled vettem.
  return 1.0 + 7.86e-4 * 101 / (celsius + 273);
}

double ndiff(double const aH) {
  return 5.6866e-7 * std::exp(0.0992063 * aH) / std::pow(0.0214465 + std::exp(0.0992063 * aH), 2.0);
}

class Sin final {
private:
  Doub eps;

public:
  Sin(Doub const aEps) : eps(aEps) {}

  void operator() (const Doub aX, VecDoub_I const &aY, VecDoub_O &aDydx) {
    aDydx[0] = (-aY[0] * ndiff(aX) / n(aX)) / eps;
  }
};

void comp(char * aPre, double const aSinIncl) {
  const Int nvar = 1;
  const Doub atol=1.0e-6, rtol=atol, h1=0.01, hmin=0.0, x1=1.0, x2=0.5;  // Minden 1 m-es magassagrol indul, es fel meterig nezem.
  VecDoub ystart(nvar);
  ystart[0]=aSinIncl;    // Kezdo beesesi szog szinusza.
  Output out(100);
  Sin sin(0.001);
  Odeint<StepperDopr853<Sin>> ode(ystart, x1, x2, atol, rtol, h1, hmin, out, sin);
  ode.integrate();
  std::cout << aPre << "_x=[";
  double xCum = 0.0;
  for (Int i=0; i<out.count; i++) {
    std::cout << xCum << (i < out.count - 1 ? ", " : "];\n");
    xCum += (out.xsave[0] - out.xsave[1]) * out.ysave[0][i] / std::sqrt(1.0 - out.ysave[0][i] * out.ysave[0][i]);  // Teljes visszaverodesnel negativ lesz, ott ez a modszer meghal.
  }
  std::cout << aPre << "h=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.xsave[i] << (i < out.count - 1 ? ", " : "];\n");
  }
  std::cout << "\n";
}

int main() {
  comp("i0999", 0.999);
  comp("i09995", 0.9995);
  comp("i09999", 0.9999);
  comp("i099995", 0.99995);
}
