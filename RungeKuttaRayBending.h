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
  struct Parameters {
    StepperType mStepper;
    double      mDistAlongRay;
    double      mTolAbs;
    double      mTolRel;
    double      mStep1;
    double      mStepMin;
    double      mStepMax;
  };

  struct Result {
    bool   mValid;
    Vertex mValue;
  };

  RungeKuttaRayBending(Parameters const &aParameters, Eikonal const &aDiffEq)
    : mDiffEq(aDiffEq)
    , mSolver(aParameters.mStepper, 0.0, aParameters.mDistAlongRay, aParameters.mTolAbs, aParameters.mTolRel,
              aParameters.mStep1, aParameters.mStepMin, aParameters.mStepMax, aDiffEq) {}

  RungeKuttaRayBending(RungeKuttaRayBending const&) = default;
  RungeKuttaRayBending(RungeKuttaRayBending &&) = delete;
  RungeKuttaRayBending& operator=(RungeKuttaRayBending const&) = delete;
  RungeKuttaRayBending& operator=(RungeKuttaRayBending &&) = delete;


  Result solve4x(Vertex const &aStart, Vector const &aDir, double const aX) {
    return mDiffEq.getEarthForm() == Eikonal::EarthForm::cFlat ? solve4xFlat(aStart, aDir, aX) : solve4xRound(aStart, aDir, aX);
  }

private:
  Result solve4xFlat(Vertex const &aStart, Vector const &aDir, double const aX);
  Result solve4xRound(Vertex const &aStart, Vector const &aDir, double const aX);
};

#endif
