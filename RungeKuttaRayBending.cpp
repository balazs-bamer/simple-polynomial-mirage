#include "RungeKuttaRayBending.h"
#include <cmath>


RungeKuttaRayBending::Result RungeKuttaRayBending::solve4xFlat(Vertex const &aStart, Vector const &aDir, double const aX) {
  typename Eikonal::Variables start;
  start[0u] = aStart(0u);
  start[1u] = aStart(1u);
  start[2u] = aStart(2u);
  auto slowness = mDiffEq.getSlowness(aStart(1u));  // from height
  start[3u] = aDir(0u) * slowness;
  start[4u] = aDir(1u) * slowness;
  start[5u] = aDir(2u) * slowness;
  auto solution = mSolver.solve(start, [aX](double const, typename Eikonal::Variables const& aY){ return aY[0] >= aX; });
  Result result;
  result.mValid = solution.mValid;
  result.mValue(0u) = solution.mValue[0u];
  result.mValue(1u) = solution.mValue[1u];
  result.mValue(2u) = solution.mValue[2u];
  return result;
}

RungeKuttaRayBending::Result RungeKuttaRayBending::solve4xRound(Vertex const &aStart, Vector const &aDir, double const aX) {
  typename Eikonal::Variables start;
  start[0u] = aStart(0u);
  start[1u] = aStart(1u) + Eikonal::csRadius;
  start[2u] = aStart(2u);
  auto slowness = mDiffEq.getSlowness(aStart(1u));  // from height
  start[3u] = aDir(0u) * slowness;
  start[4u] = aDir(1u) * slowness;
  start[5u] = aDir(2u) * slowness;
  auto solution = mSolver.solve(start, [aX](double const, typename Eikonal::Variables const& aY){ return aY[0] >= aX; });     // We now neglect the variation in perpendicular along the travelled distance.
  Result result;
  result.mValid = solution.mValid;
  result.mValue(0u) = solution.mValue[0u];
  result.mValue(1u) = solution.mValue[1u] - Eikonal::csRadius ;
  result.mValue(2u) = solution.mValue[2u];
  return result;
}


/*ShepardRayBending::Static::Static() {
  mHeightLimit = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csHeightLimit, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mB           = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csB,           std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mDelta       = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csDelta,       std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
}

ShepardRayBending::ShepardRayBending(double const aTempDiffSurface)
  : mTempDiffSurface(std::max(0.0, aTempDiffSurface))
  , mRandomDistribution(-1.0, 1.0) {
  static Static tempCoeffProfiles;

  mHeightLimit = tempCoeffProfiles.mHeightLimit->eval(csTempAmbient + mTempDiffSurface);
  mB           = tempCoeffProfiles.mB->eval(csTempAmbient + mTempDiffSurface);
  mDelta       = tempCoeffProfiles.mDelta->eval(csTempAmbient + mTempDiffSurface);
  initReflection();
}

double ShepardRayBending::getTempRiseAtHeight(double const aHeight) const {
  auto height = 100.0 * aHeight;
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - mDelta)) - csTempAmbient;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double ShepardRayBending::getHeightAtTempRise(double const aTempRise) const {
  return ::exp(::log(::log((aTempRise + csTempAmbient) / csTempAmbient) / mB) / (1.0 - mDelta)) / 100.0;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double ShepardRayBending::getRefractionAtTempRise(double const aTempRise) const {
  auto tempKelvin = csTempAmbient + aTempRise;
  auto tempCelsius = tempKelvin - csCelsius2kelvin;
  auto r = 1.0 + 7.86e-4 * csAtmosphericPressureKpa / tempKelvin - 1.5e-11 * csRelativeHumidityPercent * (tempCelsius * tempCelsius + 160.0);
  return r;
}*/
