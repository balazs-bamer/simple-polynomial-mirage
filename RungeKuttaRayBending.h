#ifndef RUNGEKUTTARAYBENDING
#define RUNGEKUTTARAYBENDING

#include "mathUtil.h"
#include "3dGeomUtil.h"
#include "Eikonal.h"
#include "OdeSolverGsl.h"


class RungeKuttaRayBending final {
private:
  Eikonal const        &mDiffEq;
  OdeSolverGsl<Eikonal> mSolver;

public:
  RungeKuttaRayBending(StepperType const aStepper, double const aDistAlongRay, double const aTolAbs, double const aTolRel, double const aStep1, Eikonal const &aDiffEq)
    : mDiffEq(aDiffEq)
    , mSolver(aStepper, 0.0, aDistAlongRay, aTolAbs, aTolRel, aStep1, aDiffEq) {}

  Vertex solve4x(Vertex const &aStart, Vector const &aDir, double const aX) {
    return mDiffEq.getEarthForm() == Eikonal::EarthForm::cFlat ? solve4xFlat(aStart, aDir, aX) : solve4xRound(aStart, aDir, aX);
  }

private:
  Vertex solve4xFlat(Vertex const &aStart, Vector const &aDir, double const aX);
  Vertex solve4xRound(Vertex const &aStart, Vector const &aDir, double const aX);
};

#endif
