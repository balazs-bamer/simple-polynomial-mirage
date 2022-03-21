#include "nr3.h"
#include "stepper.h"
#include "stepperdopr853.h"
#include "odeint.h"

class Sin final {
private:
  Doub eps;

public:
  Sin(Doub const aEps) : eps(aEps) {}

  void operator() (const Doub aX, VecDoub_I const &aY, VecDoub_O &aDydx) {
    aDydx[0] = std::cos(aX);   // dy/dx = cos(x)
  }
};

int main() {
  const Int nvar = 1;
  const Doub atol=1.0e-6, rtol=atol, h1=0.01, hmin=0.0, x1=0.0, x2=3.1415926539 * 6.0;
  VecDoub ystart(nvar);
  ystart[0]=0.0;
  Output out(-1);
  Sin sin(1.0);   // epsilon not used now
  Odeint<StepperDopr853<Sin>> ode(ystart, x1, x2, atol, rtol, h1, hmin, out, sin);
  ode.integrate();
  std::cout << "x=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.xsave[i] << ", ";
  }
  std::cout << "\ny=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.ysave[0][i] - std::sin(out.xsave[i]) << ", ";
  }
  std::cout << "\n";
}
