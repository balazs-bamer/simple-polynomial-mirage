#ifndef RUNGEKUTTARAYBENDING
#define RUNGEKUTTARAYBENDING

#include "mathUtil.h"
#include "3dGeomUtil.h"
#include "Eikonal.h"
#include "OdeSolverGsl.h"


class RungeKuttaRayBending final {
private:
  Eikonal const        &mDiffEqu;
  OdeSolverGsl<Eikonal> mSolver;

public:
  RungeKuttaRayBending(StepperType const aStepper, double const aDistAlongRay, double const aTolAbs, double const aTolRel, double const aStep1, Eikonal const &aDiffEqu)
    : mDiffEqu(aDiffEqu)
    , mSolver(aStepper, 0.0, aDistAlongRay, aTolAbs, aTolRel, aStep1, aDiffEqu) {}

  Vertex solve4x(Vertex const &aStart, Vector const &aDir, double const aX) {
    typename Eikonal::Variables start;
    start[0u] = aStart(0u);
    start[1u] = aStart(1u);
    start[2u] = aStart(2u);
    auto slowness = mDiffEqu.getSlowness(aStart(1u));  // from height
    start[3u] = aDir(0u) * slowness;
    start[4u] = aDir(1u) * slowness;
    start[5u] = aDir(2u) * slowness;
    auto [t, solution] = mSolver.solve(start, [aX](typename Eikonal::Variables const& aY){ return aY[0] >= aX; });
    Vertex result;
    result(0u) = solution[0u];
    result(1u) = solution[1u];
    result(2u) = solution[2u];
    return result;
  }
};

#endif
