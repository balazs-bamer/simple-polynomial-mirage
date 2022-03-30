#ifndef RUNGEKUTTARAYBENDING
#define RUNGEKUTTARAYBENDING

#include "mathUtil.h"
#include "3dGeomUtil.h"
#include "Eikonal.h"
#include "OdeSolver.h"
#include "StepperDormandPrice853.h"


class RungeKuttaRayBending final {
private:
  Eikonal const                             &mDiffEqu;
  OdeSolver<StepperDormandPrice853<Eikonal>> mSolver;

public:
  RungeKuttaRayBending(double const aDistAlongRay, double const aTolAbs, double const aTolRel, double const aStep1, double const aStepMin, Eikonal const &aDiffEqu)
    : mDiffEqu(aDiffEqu)
    , mSolver(0.0, aDistAlongRay, aTolAbs, aTolRel, aStep1, aStepMin, aDiffEqu) {}

  Vertex solve4x(Vertex const &aStart, Vector const &aDir, double const aX) {
    typename Eikonal::Variables start;
    start[0u] = aStart(0u);
    start[1u] = aStart(1u);
    start[2u] = aStart(2u);
    auto slowness = mDiffEqu.getSlowness(aStart(1u));  // from height
    start[3u] = aDir(0u) * slowness;
    start[4u] = aDir(1u) * slowness;
    start[5u] = aDir(2u) * slowness;
    auto solution = mSolver.solve(start, [aX](typename Eikonal::Variables const& aY){ return aY[0] >= aX; }, 0.1);
    Vertex result;
    result(0u) = solution[0u];
    result(1u) = solution[1u];
    result(2u) = solution[2u];
    return result;
  }
};

#endif
