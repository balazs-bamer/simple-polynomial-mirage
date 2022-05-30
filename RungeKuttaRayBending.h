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
  double                mMaxCosDirChange;

public:
  struct Parameters {
    StepperType mStepper;
    double      mDistAlongRay;
    double      mTolAbs;
    double      mTolRel;
    double      mStep1;
    double      mStepMin;
    double      mStepMax;
    double      mMaxCosDirChange;
  };

  struct Result {
    bool   mValid;
    Vertex mValue;
    Vector mDirection;
  };

  RungeKuttaRayBending(Parameters const &aParameters, Eikonal const &aDiffEq)
    : mDiffEq(aDiffEq)
    , mSolver(aParameters.mStepper, 0.0, aParameters.mDistAlongRay, aParameters.mTolAbs, aParameters.mTolRel,
              aParameters.mStep1, aParameters.mStepMin, aParameters.mStepMax, aDiffEq)
    , mMaxCosDirChange(aParameters.mMaxCosDirChange) {}

  RungeKuttaRayBending(RungeKuttaRayBending const&) = default;
  RungeKuttaRayBending(RungeKuttaRayBending &&) = delete;
  RungeKuttaRayBending& operator=(RungeKuttaRayBending const&) = delete;
  RungeKuttaRayBending& operator=(RungeKuttaRayBending &&) = delete;

  double getRefract(double const aH) const { return mDiffEq.getRefract(aH); }

  Result solve4x(Vertex const &aStart, Vector const &aDir, double const aX) {
    return mDiffEq.getEarthForm() == Eikonal::EarthForm::cFlat ? solve4xFlat(aStart, aDir, aX) : solve4xRound(aStart, aDir, aX);
  }

private:
  Result solve4xFlat(Vertex const &aStart, Vector const &aDir, double const aX);
  Result solve4xRound(Vertex const &aStart, Vector const &aDir, double const aX);

  bool decide2resetBigStep(typename Eikonal::Variables const& aYprev, typename Eikonal::Variables const& aYnow) {
    Vector dirPrev(aYprev[3u], aYprev[4u], aYprev[5u]);
    Vector dir(aYnow[3u], aYnow[4u], aYnow[5u]);
    auto dirChangeCos = dir.dot(dirPrev) / dir.norm() / dirPrev.norm();
    return dirChangeCos < mMaxCosDirChange;
  }
};

#endif
