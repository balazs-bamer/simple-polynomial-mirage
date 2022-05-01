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
  struct Result {
    bool   mValid;
    Vertex mValue;
  };

  RungeKuttaRayBending(StepperType const aStepper, double const aDistAlongRay, double const aTolAbs, double const aTolRel, double const aStep1, Eikonal const &aDiffEq)
    : mDiffEq(aDiffEq)
    , mSolver(aStepper, 0.0, aDistAlongRay, aTolAbs, aTolRel, aStep1, aDiffEq) {}

  Result solve4x(Vertex const &aStart, Vector const &aDir, double const aX) {
    return mDiffEq.getEarthForm() == Eikonal::EarthForm::cFlat ? solve4xFlat(aStart, aDir, aX) : solve4xRound(aStart, aDir, aX);
  }

private:
  Result solve4xFlat(Vertex const &aStart, Vector const &aDir, double const aX);
  Result solve4xRound(Vertex const &aStart, Vector const &aDir, double const aX);
};

#endif
