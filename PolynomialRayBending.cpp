#include "PolynomialRayBending.h"

#include <functional>
#include <vector>
#include <cmath>


PolynomialRayBending::Static::Static() {
  mHeightLimit = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csHeightLimit, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mB           = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csB,           std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mDelta       = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csDelta,       std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
}

PolynomialRayBending::PolynomialRayBending(double const aTempDiffSurface)
  : mTempDiffSurface(std::max(0.0, aTempDiffSurface)) {
  static Static tempCoeffProfiles;

  mHeightLimit = tempCoeffProfiles.mHeightLimit->eval(csTempAmbient + mTempDiffSurface);
  mB           = tempCoeffProfiles.mB->eval(csTempAmbient + mTempDiffSurface);
  mDelta       = tempCoeffProfiles.mDelta->eval(csTempAmbient + mTempDiffSurface);
  initReflection();
}

double PolynomialRayBending::getTempRiseAtHeight(double const aHeight) const {
  auto height = 100.0 * aHeight;
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - mDelta)) - csTempAmbient;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double PolynomialRayBending::getHeightAtTempRise(double const aTempRise) const {
  return ::exp(::log(::log((aTempRise + csTempAmbient) / csTempAmbient) / mB) / (1.0 - mDelta)) / 100.0;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double PolynomialRayBending::getRefractionAtTempRise(double const aTempRise) const {
  auto tempKelvin = csTempAmbient + aTempRise;
  auto tempCelsius = tempKelvin - csCelsius2kelvin;
  auto r = 1.0 + 7.86e-4 * csAtmosphericPressureKpa / tempKelvin - 1.5e-11 * csRelativeHumidityPercent * (tempCelsius * tempCelsius + 160.0);
  return r;
}

void PolynomialRayBending::initReflection() {
  mCriticalInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return getReflectionDispDepth(aInclination) ? -1.0 : 1.0;
  }) + csEpsilon;
  std::vector<double> inclinations;
  std::vector<double> disps;
  std::vector<double> depths;
  std::vector<double> iters;
  auto increment = (csAlmostHorizontal - mCriticalInclination) / (csInclDepthProfilePointCount - 1u);
  auto inclination = mCriticalInclination;
  for(uint32_t i = 0u; i < csInclDepthProfilePointCount; ++i) {
    inclinations.push_back(inclination);
    auto dd = getReflectionDispDepth(inclination).value();
    disps.push_back(dd.mDisp);
    depths.push_back(dd.mDepth);
    iters.push_back(dd.mIter);
    inclination += increment;
  }
  mInclination2horizDisp    = std::move(std::make_unique<PolynomApprox>(disps, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{inclinations, csInclinationProfileDegree}}));
  mInclination2virtualDepth = std::move(std::make_unique<PolynomApprox>(depths, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{inclinations, csInclinationProfileDegree}}));
  mInclination2iterations   = std::move(std::make_unique<PolynomApprox>(iters, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{inclinations, csInclinationProfileDegree}}));
}

std::optional<PolynomialRayBending::DispDepth> PolynomialRayBending::getReflectionDispDepth(double const aInclination0) const {
  double currentHeight = csInitialHeight;
  double currentTemp = getTempRiseAtHeight(currentHeight);
  double currentRefractionIndex = getRefractionAtTempRise(currentTemp);
  double sinInclination0 = ::sin(aInclination0);
  double currentSinInclination = sinInclination0;
  double horizDisp = 0.0;

  std::optional<DispDepth> result;
  uint32_t iterations = 0u;
  while(currentHeight > csMinimalHeight) {
    ++iterations;
    double nextTemp = currentTemp + static_cast<long double>(csLayerDeltaTemp);
    double nextHeight = getHeightAtTempRise(nextTemp);
    double nextRefractionIndex = getRefractionAtHeight(nextHeight);
    double factor = currentRefractionIndex / nextRefractionIndex;
    double nextSinInclination = currentSinInclination * factor;
    if(nextSinInclination >= static_cast<double>(1.0)) { // The problem here is how much sinInclanation is greater than 1. Optimal would be to have it just touch 1, so the last temperature step can't be uniform.
      auto const criticalTemp = binarySearch(currentTemp, nextTemp, csEpsilon, [this, currentRefractionIndex, currentSinInclination](double const aTemp) {
        double height = getHeightAtTempRise(aTemp);
        double refractionIndex = getRefractionAtHeight(height);
        double factor = currentRefractionIndex / refractionIndex;
        double sinInclination = currentSinInclination * factor;
        return (sinInclination > static_cast<double>(1.0) ? 1.0 : -1.0);
      });
      double criticalHeight = getHeightAtTempRise(criticalTemp);
      horizDisp += static_cast<double>(currentHeight - criticalHeight) * currentSinInclination / ::sqrt(static_cast<double>(1.0) - currentSinInclination * currentSinInclination);
      auto virt = csInitialHeight - horizDisp * ::sqrt(1.0 - sinInclination0 * sinInclination0) / sinInclination0;
      result = DispDepth(horizDisp * 2.0, virt, iterations);
      break;
    }
    else {
      horizDisp += static_cast<double>(currentHeight - nextHeight) * currentSinInclination / ::sqrt(static_cast<double>(1.0) - currentSinInclination * currentSinInclination);
      currentTemp = nextTemp;
      currentHeight = nextHeight;
      currentRefractionIndex = nextRefractionIndex;
      currentSinInclination = nextSinInclination;
    }
  }
  return result;
}
